from __future__ import print_function

import argparse
import os
import random

import chainer
import numpy as np
from chainer import training, Variable
from chainer.training import extensions

from dataset import H5pyDataset
from model import DCGAN_G, DCGAN_D, init_bn, init_conv
from sampler import sampler
from updater import WassersteinUpdater


def make_optimizer(model, lr):
    optimizer = chainer.optimizers.RMSprop(lr=lr)
    optimizer.setup(model)
    return optimizer


def main():
    parser = argparse.ArgumentParser(description='Train Unsupervised Blending GAN')
    parser.add_argument('--nz', type=int, default=100, help='Size of the latent z vector')
    parser.add_argument('--ngf', type=int, default=64, help='# of base filters in G')
    parser.add_argument('--ndf', type=int, default=64, help='# of base filters in D')
    parser.add_argument('--nc', type=int, default=3, help='# of output channels in G')
    parser.add_argument('--load_size', type=int, default=64, help='Scale image to load_size')
    parser.add_argument('--image_size', type=int, default=64, help='The height / width of the input image to network')

    parser.add_argument('--gpu', type=int, default=0, help='GPU ID (negative value indicates CPU)')
    parser.add_argument('--lr_d', type=float, default=0.00005, help='Learning rate for Critic, default=0.00005')
    parser.add_argument('--lr_g', type=float, default=0.00005, help='Learning rate for Generator, default=0.00005')
    parser.add_argument('--d_iters', type=int, default=5, help='# of D iters per each G iter')
    parser.add_argument('--n_epoch', type=int, default=25, help='# of epochs to train for')
    parser.add_argument('--clamp_lower', type=float, default=-0.01, help='Lower bound for clipping')
    parser.add_argument('--clamp_upper', type=float, default=0.01, help='Upper bound for clipping')

    parser.add_argument('--data_root', help='Path to dataset')
    parser.add_argument('--experiment', default='Wasserstein_GAN_result', help='Where to store samples and models')
    parser.add_argument('--workers', type=int, default=10, help='# of data loading workers')
    parser.add_argument('--batch_size', type=int, default=128, help='input batch size')
    parser.add_argument('--test_size', type=int, default=64, help='Batch size for testing')

    parser.add_argument('--manual_seed', type=int, default=5, help='Manul seed')

    parser.add_argument('--resume', default='', help='Resume the training from snapshot')
    parser.add_argument('--snapshot_interval', type=int, default=1, help='Interval of snapshot (epoch)')
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
    G = DCGAN_G(args.image_size, args.nc, args.ngf, init_conv, init_bn)
    print('\tInit D network ...')
    D = DCGAN_D(args.image_size, args.ndf, 1, init_conv, init_bn)
    if args.gpu >= 0:
        print('\tCopy models to gpu {} ...'.format(args.gpu))
        chainer.cuda.get_device(args.gpu).use()  # Make a specified GPU current
        G.to_gpu()  # Copy the model to the GPU
        D.to_gpu()
    print('Init models done ...\n')
    # Setup an optimizer
    optimizer_d = make_optimizer(D, args.lr_d)
    optimizer_g = make_optimizer(G, args.lr_g)

    ########################################################################################################################
    # Setup dataset & iterator
    print('Load images from {} ...'.format(args.data_root))
    trainset = H5pyDataset(args.data_root, load_size=args.load_size, crop_size=args.image_size)
    print('\tTrainset contains {} image files'.format(len(trainset)))
    print('')
    train_iter = chainer.iterators.MultiprocessIterator(trainset, args.batch_size, n_processes=args.workers,
                                                        n_prefetch=args.workers)
    ########################################################################################################################

    # Set up a trainer
    updater = WassersteinUpdater(
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
        'iteration', 'main/loss', 'D/loss', 'D/loss_real', 'D/loss_fake'
    ]), trigger=print_interval)
    trainer.extend(extensions.ProgressBar(update_interval=args.print_interval))

    trainer.extend(extensions.dump_graph('D/loss', out_name='TrainGraph.dot'))

    # Plot
    plot_interval = (args.plot_interval, 'iteration')

    trainer.extend(
        extensions.PlotReport(['main/loss'], 'iteration', file_name='loss.png', trigger=plot_interval),
        trigger=plot_interval)
    trainer.extend(
        extensions.PlotReport(['D/loss'], 'iteration', file_name='d_loss.png', trigger=plot_interval),
        trigger=plot_interval)
    trainer.extend(
        extensions.PlotReport(['D/loss_real'], 'iteration', file_name='loss_real.png', trigger=plot_interval),
        trigger=plot_interval)
    trainer.extend(
        extensions.PlotReport(['D/loss_fake'], 'iteration', file_name='loss_fake.png', trigger=plot_interval),
        trigger=plot_interval)

    # Eval
    path = os.path.join(args.experiment, 'samples')
    if not os.path.isdir(path):
        os.makedirs(path)
    print('Saving samples to {} ...\n'.format(path))

    noisev = Variable(np.asarray(np.random.normal(size=(args.test_size, args.nz, 1, 1)), dtype=np.float32))
    noisev.to_gpu(args.gpu)
    trainer.extend(sampler(G, path, noisev, 'fake_samples_{}.png'), trigger=plot_interval)

    if args.resume:
        # Resume from a snapshot
        print('Resume from {} ... \n'.format(args.resume))
        chainer.serializers.load_npz(args.resume, trainer)

    # Run the training
    print('Training start ...\n')
    trainer.run()


if __name__ == '__main__':
    main()
