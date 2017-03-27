from __future__ import print_function
import argparse

import numpy as np

from chainer import Variable, link, serializers

from chainer.links import caffe

from model import RealismCNN

def cnn2fcn(src, dst):
    print('Copying layers %s -> %s:' % (src.__class__.__name__, dst.__class__.__name__))

    for child in src.children():
        dst_child = dst['fc8' if child.name.startswith('fc8') else child.name]

        if isinstance(child, link.Link):
            print('Copying {} ...'.format(child.name))

            if child.name.startswith('fc'):
                dst_child.__dict__['W'].data[...] = np.reshape(child.__dict__['W'].data, dst_child.__dict__['W'].data.shape)
                dst_child.__dict__['b'].data[...] = np.reshape(child.__dict__['b'].data, dst_child.__dict__['b'].data.shape)
            else:
                dst_child.copyparams(child)
            
            print('\tlayer: %s -> %s' % (child.name, dst_child.name))

    return dst

def main():
    parser = argparse.ArgumentParser(description='Load caffe model for chainer')
    parser.add_argument('--caffe_model_path', default='models/realismCNN_all_iter3.caffemodel', help='Path for caffe model')
    parser.add_argument('--chainer_model_path', default='models/realismCNN_all_iter3.npz', help='Path for saving chainer model')
    args = parser.parse_args()

    print('Load caffe model from {} ...'.format(args.caffe_model_path))
    caffe_model = caffe.CaffeFunction(args.caffe_model_path)
    print('Load caffe model, DONE')

    print('\nTurn CNN into FCN, start ...\n')
    chainer_model = RealismCNN()
    chainer_model(Variable(np.zeros((1, 3, 227, 227), dtype=np.float32), volatile='on'))
    chainer_model = cnn2fcn(caffe_model, chainer_model)

    print('\nTurn CNN into FCN, DONE. Save to {} ...'.format(args.chainer_model_path))
    serializers.save_npz(args.chainer_model_path, chainer_model)


if __name__ == '__main__':
    main()