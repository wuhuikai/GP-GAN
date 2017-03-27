import os
import glob
import argparse

import numpy as np

from xml.dom import minidom

from skimage.draw import polygon

from skimage.io import imsave

def getTextByName(node, name):
    return node.getElementsByTagName(name)[0].childNodes[0].data

def main():
    parser = argparse.ArgumentParser(description='Generate mask images from xml files')
    parser.add_argument('--xml_path', help='Path for loading xml files')
    parser.add_argument('--mask_path', default='masks', help='Path for storing generated mask images')
    args = parser.parse_args()

    # Init mask folder
    if not os.path.isdir(args.mask_path):
        os.makedirs(args.mask_path)
    print('Masks will save to {} ...\n'.format(args.mask_path))

    xmls = glob.glob(os.path.join(args.xml_path, '*.xml'))
    total_size = len(xmls)

    for idx, xml_path in enumerate(xmls):
        name = os.path.splitext(os.path.basename(xml_path))[0]
        print('Processing {}/{}, name = {} ...'.format(idx+1, total_size, name))
        if not os.path.isdir(os.path.join(args.mask_path, name)):
        	os.makedirs(os.path.join(args.mask_path, name))
        print('\tSave to subfloder {} ...'.format(name))

        dom = minidom.parse(open(xml_path))
        image_size_ele = dom.getElementsByTagName('imagesize')[0]
        size = (int(getTextByName(image_size_ele, 'nrows')), int(getTextByName(image_size_ele, 'ncols')))

        i = 0
        for obj in dom.getElementsByTagName('object'):
            obj_name = getTextByName(obj, 'name')
            if int(getTextByName(obj, 'deleted')):
            	continue
            
            polygon_ele = obj.getElementsByTagName('polygon')[0]
            pts_ele = polygon_ele.getElementsByTagName('pt')
            r = np.array([int(getTextByName(pt, 'y')) for pt in pts_ele])
            c = np.array([int(getTextByName(pt, 'x')) for pt in pts_ele])
            rr, cc = polygon(r, c)

            mask = np.zeros(size, dtype=np.uint8)
            mask[rr, cc] = 1
            
            imsave(os.path.join(args.mask_path, name, '{}_{}.png'.format(obj_name, i)), mask)
            i += 1

if __name__ == '__main__':
    main()