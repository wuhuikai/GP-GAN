from __future__ import print_function

import argparse
import os
import random

import matplotlib

matplotlib.use('Agg')

import chainer
from chainer import training, Variable
from chainer.training import extensions

from model import EncoderDecoder, DCGAN_D, init_bn, init_conv
from dataset import BlendingDataset
from updater import EncoderDecoderBlendingUpdater
from sampler import sampler


def make_optimizer(model, alpha, beta1):
    optimizer = chainer.optimizers.Adam(alpha=alpha, beta1=beta1)
    optimizer.setup(model)
    return optimizer


def main():
    parser = argparse.ArgumentParser(description='Train Blending GAN')
    parser.add_argument('--nef', type=int, default=64, help='# of base filters in encoder')
    parser.add_argument('--ngf', type=int, default=64, help='# of base filters in decoder')
    parser.add_argument('--nc', type=int, default=3, help='# of output channels in decoder')
    parser.add_argument('--nBottleneck', type=int, default=4000, help='# of output channels in encoder')
    parser.add_argument('--ndf', type=int, default=64, help='# of base filters in D')

    parser.add_argument('--lr_d', type=float, default=0.0002, help='Learning rate for Critic, default=0.0002')
    parser.add_argument('--lr_g', type=float, default=0.002, help='Learning rate for Generator, default=0.002')
    parser.add_argument('--beta1', type=float, default=0.5, help='Beta for Adam, default=0.5')
    parser.add_argument('--l2_weight', type=float, default=0.999, help='Weight for l2 loss, default=0.999')

    parser.add_argument('--gpu', type=int, default=0, help='GPU ID (negative value indicates CPU)')
    parser.add_argument('--n_epoch', type=int, default=25, help='# of epochs to train for')

    parser.add_argument('--data_root', help='Path to dataset')
    parser.add_argument('--load_size', type=int, default=64, help='Scale image to load_size')
    parser.add_argument('--image_size', type=int, default=64, help='The height / width of the input image to network')
    parser.add_argument('--ratio', type=float, default=0.5, help='Ratio for center square size v.s. image_size')
    parser.add_argument('--val_ratio', type=float, default=0.05, help='Ratio for validation set v.s. data set')

    parser.add_argument('--d_iters', type=int, default=5, help='# of D iters per each G iter')
    parser.add_argument('--clamp_lower', type=float, default=-0.01, help='Lower bound for clipping')
    parser.add_argument('--clamp_upper', type=float, default=0.01, help='Upper bound for clipping')

    parser.add_argument('--experiment', default='encoder_decoder_blending_result',
                        help='Where to store samples and models')
    parser.add_argument('--test_folder', default='samples', help='Where to store test results')
    parser.add_argument('--workers', type=int, default=4, help='# of data loading workers')
    parser.add_argument('--batch_size', type=int, default=64, help='Input batch size')
    parser.add_argument('--test_size', type=int, default=64, help='Batch size for testing')

    parser.add_argument('--train_samples', type=int, default=150000, help='# of training examples')
    parser.add_argument('--test_samples', type=int, default=256, help='# of testing examples')

    parser.add_argument('--manual_seed', type=int, default=5, help='Manul seed')

    parser.add_argument('--resume', default='', help='Resume the training from snapshot')
    parser.add_argument('--snapshot_interval', type=int, default=1, help='Interval of snapshot (epochs)')
    parser.add_argument('--print_interval', type=int, default=1, help='Interval of printing log to console (iteration)')
    parser.add_argument('--plot_interval', type=int, default=10, help='Interval of plot (iteration)')
    args = parser.parse_args()

    random.seed(args.manual_seed)

    print('Input arguments:')
    for key, value in vars(args).items():
        print('\t{}: {}'.format(key, value))
    print('')

    # Set up G & D
    print('Create & Init models ...')
    print('\tInit G network ...')
    G = EncoderDecoder(args.nef, args.ngf, args.nc, args.nBottleneck, image_size=args.image_size, conv_init=init_conv,
                       bn_init=init_bn)
    print('\tInit D network ...')
    D = DCGAN_D(args.image_size, args.ndf, conv_init=init_conv, bn_init=init_bn)
    if args.gpu >= 0:
        print('\tCopy models to gpu {} ...'.format(args.gpu))
        chainer.cuda.get_device(args.gpu).use()  # Make a specified GPU current
        G.to_gpu()  # Copy the model to the GPU
        D.to_gpu()
    print('Init models done ...\n')
    # Setup an optimizer
    optimizer_d = make_optimizer(D, args.lr_d, args.beta1)
    optimizer_g = make_optimizer(G, args.lr_g, args.beta1)

    ########################################################################################################################
    # Setup dataset & iterator
    print('Load images from {} ...'.format(args.data_root))
    folders = sorted(
        [folder for folder in os.listdir(args.data_root) if os.path.isdir(os.path.join(args.data_root, folder))])
    val_end = int(args.val_ratio * len(folders))
    print('\t{} folders in total, {} val folders ...'.format(len(folders), val_end))
    trainset = BlendingDataset(args.train_samples, folders[val_end:], args.data_root, args.ratio, args.load_size,
                               args.image_size)
    valset = BlendingDataset(args.test_samples, folders[:val_end], args.data_root, args.ratio, args.load_size,
                             args.image_size)
    print('\tTrainset contains {} image files'.format(len(trainset)))
    print('\tValset contains {} image files'.format(len(valset)))
    print('')
    train_iter = chainer.iterators.MultiprocessIterator(trainset, args.batch_size, n_processes=args.workers,
                                                        n_prefetch=args.workers)
    ########################################################################################################################

    # Set up a trainer
    updater = EncoderDecoderBlendingUpdater(
        models=(G, D),
        args=args,
        iterator=train_iter,
        optimizer={'main': optimizer_g, 'D': optimizer_d},
        device=args.gpu
    )
    trainer = training.Trainer(updater, (args.n_epoch, 'epoch'), out=args.experiment)

    # Snapshot
    snapshot_interval = (args.snapshot_interval, 'epoch')
    trainer.extend(
        extensions.snapshot(filename='snapshot_epoch_{.updater.epoch}.npz'),
        trigger=snapshot_interval)
    trainer.extend(extensions.snapshot_object(
        G, 'g_epoch_{.updater.epoch}.npz'), trigger=snapshot_interval)
    trainer.extend(extensions.snapshot_object(
        D, 'd_epoch_{.updater.epoch}.npz'), trigger=snapshot_interval)

    # Display
    print_interval = (args.print_interval, 'iteration')
    trainer.extend(extensions.LogReport(trigger=print_interval))
    trainer.extend(extensions.PrintReport([
        'iteration', 'main/loss', 'D/loss', 'main/l2_loss'
    ]), trigger=print_interval)
    trainer.extend(extensions.ProgressBar(update_interval=args.print_interval))

    trainer.extend(extensions.dump_graph('D/loss', out_name='TrainGraph.dot'))

    # Plot
    plot_interval = (args.plot_interval, 'iteration')

    trainer.extend(
        extensions.PlotReport(['main/loss'], 'iteration', file_name='loss.png', trigger=plot_interval))
    trainer.extend(
        extensions.PlotReport(['D/loss'], 'iteration', file_name='d_loss.png', trigger=plot_interval))
    trainer.extend(
        extensions.PlotReport(['main/l2_loss'], 'iteration', file_name='l2_loss.png', trigger=plot_interval))

    # Eval
    path = os.path.join(args.experiment, args.test_folder)
    if not os.path.isdir(path):
        os.makedirs(path)
    print('Saving samples to {} ...\n'.format(path))

    train_batch = [trainset[idx][0] for idx in range(args.test_size)]
    train_v = Variable(chainer.dataset.concat_examples(train_batch, args.gpu))
    trainer.extend(sampler(G, path, train_v, 'fake_samples_train_{}.png'), trigger=plot_interval)

    val_batch = [valset[idx][0] for idx in range(args.test_size)]
    val_v = Variable(chainer.dataset.concat_examples(val_batch, args.gpu))
    trainer.extend(sampler(G, path, val_v, 'fake_samples_val_{}.png'), trigger=plot_interval)

    if args.resume:
        # Resume from a snapshot
        print('Resume from {} ... \n'.format(args.resume))
        chainer.serializers.load_npz(args.resume, trainer)

    # Run the training
    print('Training start ...\n')
    trainer.run()


if __name__ == '__main__':
    main()
