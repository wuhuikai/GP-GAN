import os
import glob
import argparse

import numpy as np

from skimage.io import imread, imsave

def main():
    parser = argparse.ArgumentParser(description='Mask cropping')
    parser.add_argument('--mask_path', help='Path for masks')
    parser.add_argument('--crop_path', default='cropped_masks', help='Path for storing cropped masks')
    parser.add_argument('--bbox_path', default='../DataBase/TransientAttributes/bbox.txt', help='Path for bounding-box txt')
    args = parser.parse_args()

    # Init mask folder
    if not os.path.isdir(args.crop_path):
        os.makedirs(args.crop_path)
    print('Masks will save to {} ...\n'.format(args.crop_path))

    with open(args.bbox_path) as f:
        for line in f:
            name, bbox = line.strip().split(':')
            sx, sy, ex, ey = [int(i) for i in bbox.split(',')]

            print('Processing {} ...'.format(name))
            masks = glob.glob(os.path.join(args.mask_path, name, '*.png'))
            if not os.path.isdir(os.path.join(args.crop_path, name)):
                os.makedirs(os.path.join(args.crop_path, name))

            for mask_path in masks:
                mask = imread(mask_path)
                cropped_mask = mask[sx:ex, sy:ey]
                
                mask_name = os.path.basename(mask_path)
                imsave(os.path.join(args.crop_path, name, mask_name), cropped_mask)

if __name__ == '__main__':
    main()