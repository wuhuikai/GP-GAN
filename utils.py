import math

import numpy as np


def make_grid(tensor, padding=2):
    """
    Given a 4D mini-batch Tensor of shape (B x C x H x W), makes a grid of images
    """
    # make the mini-batch of images into a grid
    nmaps = tensor.shape[0]
    xmaps = int(nmaps ** 0.5)
    ymaps = int(math.ceil(nmaps / xmaps))
    height, width = int(tensor.shape[2] + padding), int(tensor.shape[3] + padding)
    grid = np.ones((3, height * ymaps, width * xmaps))
    k = 0
    sy = 1 + padding // 2
    for y in range(ymaps):
        sx = 1 + padding // 2
        for x in range(xmaps):
            if k >= nmaps:
                break
            grid[:, sy:sy + height - padding, sx:sx + width - padding] = tensor[k]
            sx += width
            k = k + 1
        sy += height
    return grid
