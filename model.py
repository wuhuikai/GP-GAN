import chainer
import chainer.links as L
import chainer.functions as F

from chainer import cuda

def init_conv(array):
    xp = cuda.get_array_module(array)
    array[...] = xp.random.normal(loc=0.0, scale=0.02, size=array.shape)
def init_bn(array):
    xp = cuda.get_array_module(array)
    array[...] = xp.random.normal(loc=1.0, scale=0.02, size=array.shape)

class ReLU(chainer.Chain):
    def __init__(self):
        super(ReLU, self).__init__()

    def __call__(self, x):
        return F.relu(x)

class Tanh(chainer.Chain):
    def __init__(self):
        super(Tanh, self).__init__()

    def __call__(self, x):
        return F.tanh(x)

class LeakyReLU(chainer.Chain):
    def __init__(self):
        super(LeakyReLU, self).__init__()

    def __call__(self, x):
        return F.leaky_relu(x)

class DCGAN_G(chainer.ChainList):
    def __init__(self, isize, nc, ngf, conv_init=None, bn_init=None):
        cngf, tisize = ngf//2, 4
        while tisize != isize:
            cngf = cngf * 2
            tisize = tisize * 2

        layers = []
        # input is Z, going into a convolution
        layers.append(L.Deconvolution2D(None, cngf, ksize=4, stride=1, pad=0, initialW=conv_init, nobias=True))
        layers.append(L.BatchNormalization(cngf, initial_gamma=bn_init))
        layers.append(ReLU())
        csize, cndf = 4, cngf
        while csize < isize//2:
            layers.append(L.Deconvolution2D(None, cngf//2, ksize=4, stride=2, pad=1, initialW=conv_init, nobias=True))
            layers.append(L.BatchNormalization(cngf//2, initial_gamma=bn_init))
            layers.append(ReLU())
            cngf = cngf // 2
            csize = csize * 2
        layers.append(L.Deconvolution2D(None, nc, ksize=4, stride=2, pad=1, initialW=conv_init, nobias=True))
        layers.append(Tanh())

        super(DCGAN_G, self).__init__(*layers)

    def __call__(self, x):
        for i in range(len(self)):
            x = self[i](x)

        return x

class DCGAN_D(chainer.ChainList):
    def __init__(self, isize, ndf, nz=1, conv_init=None, bn_init=None):
        layers = []
        layers.append(L.Convolution2D(None, ndf, ksize=4, stride=2, pad=1, initialW=conv_init, nobias=True))
        layers.append(LeakyReLU())
        csize, cndf = isize / 2, ndf
        while csize > 4:
            in_feat = cndf
            out_feat = cndf * 2
            layers.append(L.Convolution2D(None, out_feat, ksize=4, stride=2, pad=1, initialW=conv_init, nobias=True))
            layers.append(L.BatchNormalization(out_feat, initial_gamma=bn_init))
            layers.append(LeakyReLU())

            cndf = cndf * 2
            csize = csize / 2
        # state size. K x 4 x 4
        layers.append(L.Convolution2D(None, nz, ksize=4, stride=1, pad=0, initialW=conv_init, nobias=True))

        super(DCGAN_D, self).__init__(*layers)

    def encode(self, x):
        for i in range(len(self)):
            x = self[i](x)

        return x

    def __call__(self, x):
        x  = self.encode(x)
        x = F.sum(x, axis=0) / x.shape[0]
        return F.squeeze(x)

class EncoderDecoder(chainer.Chain):
    def __init__(self, nef, ngf, nc, nBottleneck, image_size=64, conv_init=None, bn_init=None):
        super(EncoderDecoder, self).__init__(
            encoder = DCGAN_D(image_size, nef, nBottleneck, conv_init, bn_init),
            bn      = L.BatchNormalization(nBottleneck, initial_gamma=bn_init),
            decoder = DCGAN_G(image_size, nc, ngf, conv_init, bn_init)
        )

    def encode(self, x):
        h = self.encoder.encode(x)
        h = F.leaky_relu(self.bn(h))

        return h

    def decode(self, x):
        h = self.decoder(x)

        return h

    def __call__(self, x):
        h = self.encode(x)
        h = self.decode(h)
        return h

class RealismCNN(chainer.Chain):
    def __init__(self, w_init=None):
        super(RealismCNN, self).__init__(
            conv1_1=L.Convolution2D(None, 64, ksize=3, stride=1, pad=1, initialW=w_init),
            conv1_2=L.Convolution2D(None, 64, ksize=3, stride=1, pad=1, initialW=w_init),

            conv2_1=L.Convolution2D(None, 128, ksize=3, stride=1, pad=1, initialW=w_init),
            conv2_2=L.Convolution2D(None, 128, ksize=3, stride=1, pad=1, initialW=w_init),

            conv3_1=L.Convolution2D(None, 256, ksize=3, stride=1, pad=1, initialW=w_init),
            conv3_2=L.Convolution2D(None, 256, ksize=3, stride=1, pad=1, initialW=w_init),
            conv3_3=L.Convolution2D(None, 256, ksize=3, stride=1, pad=1, initialW=w_init),

            conv4_1=L.Convolution2D(None, 512, ksize=3, stride=1, pad=1, initialW=w_init),
            conv4_2=L.Convolution2D(None, 512, ksize=3, stride=1, pad=1, initialW=w_init),
            conv4_3=L.Convolution2D(None, 512, ksize=3, stride=1, pad=1, initialW=w_init),

            conv5_1=L.Convolution2D(None, 512, ksize=3, stride=1, pad=1, initialW=w_init),
            conv5_2=L.Convolution2D(None, 512, ksize=3, stride=1, pad=1, initialW=w_init),
            conv5_3=L.Convolution2D(None, 512, ksize=3, stride=1, pad=1, initialW=w_init),

            fc6=L.Convolution2D(None, 4096, ksize=7, stride=1, pad=0, initialW=w_init),
            fc7=L.Convolution2D(None, 4096, ksize=1, stride=1, pad=0, initialW=w_init),
            fc8=L.Convolution2D(None, 2, ksize=1, stride=1, pad=0, initialW=w_init)
        )

    def __call__(self, x, dropout=True):
        h = F.relu(self.conv1_1(x))
        h = F.relu(self.conv1_2(h))
        h = F.max_pooling_2d(h, ksize=2, stride=2)

        h = F.relu(self.conv2_1(h))
        h = F.relu(self.conv2_2(h))
        h = F.max_pooling_2d(h, ksize=2, stride=2)

        h = F.relu(self.conv3_1(h))
        h = F.relu(self.conv3_2(h))
        h = F.relu(self.conv3_3(h))
        h = F.max_pooling_2d(h, ksize=2, stride=2)

        h = F.relu(self.conv4_1(h))
        h = F.relu(self.conv4_2(h))
        h = F.relu(self.conv4_3(h))
        h = F.max_pooling_2d(h, ksize=2, stride=2)

        h = F.relu(self.conv5_1(h))
        h = F.relu(self.conv5_2(h))
        h = F.relu(self.conv5_3(h))
        h = F.max_pooling_2d(h, ksize=2, stride=2)

        h = F.dropout(F.relu(self.fc6(h)), train=dropout)
        h = F.dropout(F.relu(self.fc7(h)), train=dropout)
        h = self.fc8(h)

        return h