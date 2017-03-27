from __future__ import print_function
import argparse

import chainer
from chainer import Variable, serializers

from model import RealismCNN

from utils import im_preprocess_vgg

import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Predict a list of images wheather realistic or not')
    parser.add_argument('--gpu', type=int, default=0, help='GPU ID (negative value indicates CPU)')
    parser.add_argument('--model_path', default='models/realismCNN_all_iter3.npz', help='Path for pretrained model')
    parser.add_argument('--list_path', help='Path for file storing image list')
    parser.add_argument('--batch_size', type=int, default=10, help='Batchsize of 1 iteration')
    parser.add_argument('--load_size', type=int, default=256, help='Scale image to load_size')
    parser.add_argument('--result_path', default='result.txt', help='Path for file storing results')
    args = parser.parse_args()

    model = RealismCNN()
    print('Load pretrained model from {} ...'.format(args.model_path))
    serializers.load_npz(args.model_path, model)
    if args.gpu >= 0:
        chainer.cuda.get_device(args.gpu).use()  # Make a specified GPU current
        model.to_gpu()                           # Copy the model to the GPU

    print('Load images from {} ...'.format(args.list_path))
    dataset = chainer.datasets.ImageDataset(paths=args.list_path, root='')
    print('{} images in total loaded'.format(len(dataset)))
    data_iterator = chainer.iterators.SerialIterator(dataset, args.batch_size, repeat=False, shuffle=False)

    scores = np.zeros((0, 2))
    for idx, batch in enumerate(data_iterator):
        print('Processing batch {}->{}/{} ...'.format(idx*args.batch_size+1, min(len(dataset), (idx+1)*args.batch_size), len(dataset)))
        batch = [im_preprocess_vgg(np.transpose(im, [1, 2, 0]), args.load_size) for im in batch]
        batch = Variable(chainer.dataset.concat_examples(batch, args.gpu), volatile='on')
        result = chainer.cuda.to_cpu(model(batch, dropout=False).data)
        scores = np.vstack((scores, np.mean(result, axis=(2, 3))))

    print('Processing DONE !')
    print('Saving result to {} ...'.format(args.result_path))
    with open(args.result_path, 'w') as f:
        for score in scores:
            f.write('{},{}\n'.format(score[0], score[1]))

if __name__ == '__main__':
    main()