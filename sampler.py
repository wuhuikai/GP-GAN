import os

import chainer
import numpy as np
from skimage.io import imsave

from utils import make_grid


def sampler(G, dst, inputv, name):
    @chainer.training.make_extension()
    def make_image(trainer):
        with chainer.using_config("Train", False):
            with chainer.no_backprop_mode():
                fake = G(inputv)
                fake = chainer.cuda.to_cpu(fake.data)
                img = make_grid(fake)
                img = np.asarray(np.transpose(np.clip((img + 1) * 127.5, 0, 255), (1, 2, 0)), dtype=np.uint8)
                imsave(os.path.join(dst, name.format(trainer.updater.iteration)), img)

    return make_image
