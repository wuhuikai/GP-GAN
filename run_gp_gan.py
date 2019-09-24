import argparse
import os

import chainer
from chainer import cuda, serializers
from skimage import img_as_float
from skimage.io import imread, imsave

from gp_gan import gp_gan
from model import EncoderDecoder, DCGAN_G

basename = lambda path: os.path.splitext(os.path.basename(path))[0]

"""
    Note: source image, destination image and mask image have the same size.
"""


def main():
    parser = argparse.ArgumentParser(description='Gaussian-Poisson GAN for high-resolution image blending')
    parser.add_argument('--nef', type=int, default=64, help='# of base filters in encoder')
    parser.add_argument('--ngf', type=int, default=64, help='# of base filters in decoder or G')
    parser.add_argument('--nc', type=int, default=3, help='# of output channels in decoder or G')
    parser.add_argument('--nBottleneck', type=int, default=4000, help='# of output channels in encoder')
    parser.add_argument('--ndf', type=int, default=64, help='# of base filters in D')

    parser.add_argument('--image_size', type=int, default=64, help='The height / width of the input image to network')

    parser.add_argument('--color_weight', type=float, default=1, help='Color weight')
    parser.add_argument('--sigma', type=float, default=0.5,
                        help='Sigma for gaussian smooth of Gaussian-Poisson Equation')
    parser.add_argument('--gradient_kernel', type=str, default='normal', help='Kernel type for calc gradient')
    parser.add_argument('--smooth_sigma', type=float, default=1, help='Sigma for gaussian smooth of Laplacian pyramid')

    parser.add_argument('--supervised', type=lambda x: x == 'True', default=True,
                        help='Use unsupervised Blending GAN if False')
    parser.add_argument('--nz', type=int, default=100, help='Size of the latent z vector')
    parser.add_argument('--n_iteration', type=int, default=1000, help='# of iterations for optimizing z')

    parser.add_argument('--gpu', type=int, default=0, help='GPU ID (negative value indicates CPU)')
    parser.add_argument('--g_path', default='models/blending_gan.npz', help='Path for pretrained Blending GAN model')
    parser.add_argument('--unsupervised_path', default='models/unsupervised_blending_gan.npz',
                        help='Path for pretrained unsupervised Blending GAN model')
    parser.add_argument('--list_path', default='',
                        help='File for input list in csv format: obj_path;bg_path;mask_path in each line')
    parser.add_argument('--result_folder', default='blending_result', help='Name for folder storing results')

    parser.add_argument('--src_image', default='', help='Path for source image')
    parser.add_argument('--dst_image', default='', help='Path for destination image')
    parser.add_argument('--mask_image', default='', help='Path for mask image')
    parser.add_argument('--blended_image', default='', help='Where to save blended image')

    args = parser.parse_args()

    print('Input arguments:')
    for key, value in vars(args).items():
        print('\t{}: {}'.format(key, value))
    print('')

    # Init CNN model
    if args.supervised:
        G = EncoderDecoder(args.nef, args.ngf, args.nc, args.nBottleneck, image_size=args.image_size)
        print('Load pretrained Blending GAN model from {} ...'.format(args.g_path))
        serializers.load_npz(args.g_path, G)
    else:
        chainer.config.use_cudnn = 'never'
        G = DCGAN_G(args.image_size, args.nc, args.ngf)
        print('Load pretrained unsupervised Blending GAN model from {} ...'.format(args.unsupervised_path))
        serializers.load_npz(args.unsupervised_path, G)

    if args.gpu >= 0:
        cuda.get_device(args.gpu).use()  # Make a specified GPU current
        G.to_gpu()  # Copy the model to the GPU

    # Init image list
    if args.list_path:
        print('Load images from {} ...'.format(args.list_path))
        with open(args.list_path) as f:
            test_list = [line.strip().split(';') for line in f]
        print('\t {} images in total ...\n'.format(len(test_list)))
    else:
        test_list = [(args.src_image, args.dst_image, args.mask_image)]

    if not args.blended_image:
        # Init result folder
        if not os.path.isdir(args.result_folder):
            os.makedirs(args.result_folder)
        print('Result will save to {} ...\n'.format(args.result_folder))

    total_size = len(test_list)
    for idx in range(total_size):
        print('Processing {}/{} ...'.format(idx + 1, total_size))

        # load image
        obj = img_as_float(imread(test_list[idx][0]))
        bg = img_as_float(imread(test_list[idx][1]))
        mask = imread(test_list[idx][2], as_gray=True).astype(obj.dtype)

        with chainer.using_config("train", False):
            blended_im = gp_gan(obj, bg, mask, G, args.image_size, args.gpu, color_weight=args.color_weight,
                                sigma=args.sigma,
                                gradient_kernel=args.gradient_kernel, smooth_sigma=args.smooth_sigma,
                                supervised=args.supervised,
                                nz=args.nz, n_iteration=args.n_iteration)

        if args.blended_image:
            imsave(args.blended_image, blended_im)
        else:
            imsave('{}/obj_{}_bg_{}_mask_{}.png'.format(args.result_folder, basename(test_list[idx][0]),
                                                        basename(test_list[idx][1]), basename(test_list[idx][2])),
                   blended_im)


if __name__ == '__main__':
    main()
