import argparse
import glob
import os

from skimage.io import imread, imsave


def main():
    parser = argparse.ArgumentParser(description='Cropping aligned images')
    parser.add_argument('--data_root', help='Path for aligned images')
    parser.add_argument('--result_folder', default='cropped_images', help='Path for cropped images')
    parser.add_argument('--bbox_path', default='DataBase/TransientAttributes/bbox.txt',
                        help='Path for bounding-box txt')
    args = parser.parse_args()

    # Init mask folder
    if not os.path.isdir(args.result_folder):
        os.makedirs(args.result_folder)
    print('Cropped images will be saved to {} ...\n'.format(args.result_folder))

    with open(args.bbox_path) as f:
        for line in f:
            name, bbox = line.strip().split(':')
            sx, sy, ex, ey = [int(i) for i in bbox.split(',')]

            print('Processing {} ...'.format(name))
            images = glob.glob(os.path.join(args.data_root, name, '*'))
            if not os.path.isdir(os.path.join(args.result_folder, name)):
                os.makedirs(os.path.join(args.result_folder, name))

            for data_root in images:
                mask = imread(data_root)
                cropped_mask = mask[sx:ex, sy:ey]

                mask_name = os.path.basename(data_root)
                imsave(os.path.join(args.result_folder, name, mask_name), cropped_mask)


if __name__ == '__main__':
    main()
