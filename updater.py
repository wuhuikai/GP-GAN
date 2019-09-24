import chainer
import chainer.functions as F
import numpy as np
from chainer import Variable


class WassersteinUpdaterFramework(chainer.training.StandardUpdater):
    def __init__(self, *args, **kwargs):
        self.G, self.D = kwargs.pop('models')
        self.args = kwargs.pop('args')

        super(WassersteinUpdaterFramework, self).__init__(*args, **kwargs)

    def d_loss(self, errD_real, errD_fake):
        errD = errD_real - errD_fake

        chainer.report({'loss_real': errD_real}, self.D)
        chainer.report({'loss_fake': errD_fake}, self.D)
        chainer.report({'loss': errD}, self.D)

        return errD

    def update_d(self, optimizer):
        raise NotImplementedError

    def update_g(self, optimizer):
        raise NotImplementedError

    def update_core(self):
        d_optimizer = self.get_optimizer('D')
        g_optimizer = self.get_optimizer('main')
        ############################
        # (1) Update D network
        ###########################
        # train the discriminator Diters times
        if self.iteration < 25 or self.iteration % 500 == 0:
            Diters = 100
        else:
            Diters = self.args.d_iters

        for _ in range(Diters):
            # clamp parameters to a cube
            for p in self.D.params():
                if p.data is None:
                    continue
                p.data.clip(self.args.clamp_lower, self.args.clamp_upper, p.data)

            self.update_d(d_optimizer)

        ############################
        # (2) Update G network
        ###########################
        self.update_g(g_optimizer)


class EncoderDecoderBlendingUpdater(WassersteinUpdaterFramework):
    def __init__(self, *args, **kwargs):
        super(EncoderDecoderBlendingUpdater, self).__init__(*args, **kwargs)

    def g_loss(self, errG, fake, gtv):
        l2_loss = F.mean_squared_error(fake, gtv)
        loss = (1 - self.args.l2_weight) * errG + self.args.l2_weight * l2_loss

        chainer.report({'loss': loss}, self.G)
        chainer.report({'l2_loss': l2_loss}, self.G)
        chainer.report({'gan_loss': errG}, self.G)

        return loss

    def update_d(self, optimizer):
        batch = self.get_iterator('main').next()
        inputv = Variable(self.converter([inputs for inputs, _ in batch], self.device))
        gtv = Variable(self.converter([gt for _, gt in batch], self.device))
        errD_real = self.D(gtv)

        # train with fake
        fake = self.G(inputv)
        errD_fake = self.D(fake)

        optimizer.update(self.d_loss, errD_real, errD_fake)

    def update_g(self, optimizer):
        batch = self.get_iterator('main').next()
        inputv = Variable(self.converter([inputs for inputs, _ in batch], self.device))
        gtv = Variable(self.converter([gt for _, gt in batch], self.device))
        fake = self.G(inputv)
        errG = self.D(fake)
        optimizer.update(self.g_loss, errG, fake, gtv)


class WassersteinUpdater(WassersteinUpdaterFramework):
    def __init__(self, *args, **kwargs):
        super(WassersteinUpdater, self).__init__(*args, **kwargs)

    def g_loss(self, errG):
        chainer.report({'loss': errG}, self.G)

        return errG

    def update_d(self, optimizer):
        batch = self.get_iterator('main').next()
        inputv = Variable(self.converter(batch, self.device))
        errD_real = self.D(inputv)

        # train with fake
        noisev = Variable(
            np.asarray(np.random.normal(size=(self.args.batch_size, self.args.nz, 1, 1)), dtype=np.float32))
        noisev.to_device(self.device)
        fake = self.G(noisev)
        errD_fake = self.D(fake)

        optimizer.update(self.d_loss, errD_real, errD_fake)

    def update_g(self, optimizer):
        noisev = Variable(
            np.asarray(np.random.normal(size=(self.args.batch_size, self.args.nz, 1, 1)), dtype=np.float32))
        noisev.to_device(self.device)
        fake = self.G(noisev)
        errG = self.D(fake)
        optimizer.update(self.g_loss, errG)
