# GP-GAN: Towards Realistic High-Resolution Image Blending (ACMMM 2019, **oral**)
[[Project]](https://wuhuikai.github.io/GP-GAN-Project/)   [[Paper]](https://arxiv.org/abs/1703.07195) [[Demo]](http://wuhuikai.me/DeepJS/)   [[Related Work: A2RL (for Auto Image Cropping)]](https://github.com/wuhuikai/TF-A2RL)  [[Colab]](https://colab.research.google.com/drive/11ESQZpS-Abq87WZSfc325QwCFD-j_VnL)  
Official Chainer implementation of GP-GAN: Towards Realistic High-Resolution Image Blending

## Overview

| source | destination | mask | composited | blended |
| --- | --- | --- | --- | --- |
| ![](images/test_images/src.jpg) | ![](images/test_images/dst.jpg) | ![](images/test_images/mask_display.png) | ![](images/test_images/copy-paste.png) | ![](images/test_images/result.png) |

The author's implementation of GP-GAN, the high-resolution image blending algorithm described in:  
"GP-GAN: Towards Realistic High-Resolution Image Blending"   
Huikai Wu, Shuai Zheng, Junge Zhang, Kaiqi Huang

Given a mask, our algorithm can blend the source image and the destination image, generating a high-resolution and realsitic blended image. Our algorithm is based on deep generative models [Wasserstein GAN](https://arxiv.org/abs/1701.07875).

Contact: Hui-Kai Wu (huikaiwu@icloud.com)

## Citation
```
@article{wu2017gp,
  title   = {GP-GAN: Towards Realistic High-Resolution Image Blending},
  author  = {Wu, Huikai and Zheng, Shuai and Zhang, Junge and Huang, Kaiqi},
  journal = {ACMMM},
  year    = {2019}
}
```

## Getting started
* The code is tested with `python==3.5` and `chainer==6.3.0` on `Ubuntu 16.04 LTS`.
* Download the code from GitHub:
    ```bash
    git clone https://github.com/wuhuikai/GP-GAN.git
    cd GP-GAN
    ```
* Install the requirements:
    ```bash
    pip install -r requirements/test/requirements.txt
    ```
* Download the pretrained model `blending_gan.npz` or `unsupervised_blending_gan.npz` from [Google Drive](https://drive.google.com/open?id=0Bybnpq8dvwudVjBHNWNHUmVSV28), and then put them in the folder `models`.

* Run the script for `blending_gan.npz`:
    ``` bash
    python run_gp_gan.py --src_image images/test_images/src.jpg --dst_image images/test_images/dst.jpg --mask_image images/test_images/mask.png --blended_image images/test_images/result.png
    ```
    **Or** run the script for `unsupervised_blending_gan.npz`:
    ``` bash
    python run_gp_gan.py --src_image images/test_images/src.jpg --dst_image images/test_images/dst.jpg --mask_image images/test_images/mask.png --blended_image images/test_images/result.png --supervised False
    ```
* Type `python run_gp_gan.py --help` for a complete list of the arguments.

## Train GP-GAN step by step
### Train Blending GAN
* Download Transient Attributes Dataset [here](http://transattr.cs.brown.edu/files/aligned_images.tar).
* Crop the images in each subfolder:
    ```bash
    python crop_aligned_images.py --data_root [Path for imageAlignedLD in Transient Attributes Dataset]
    ```
* Train Blending GAN:
    ```bash
    python train_blending_gan.py --data_root [Path for cropped aligned images of Transient Attributes Dataset]
    ```
* Training Curve

    ![](images/blending_gan_result/loss.png)
* Visual Result

    | Training Set | Validation Set |
    | --- | --- |
    | ![](images/blending_gan_result/train.png) | ![](images/blending_gan_result/val.png) |

### Training Unsupervised Blending GAN
* Requirements
    ```bash
    pip install git+git://github.com/mila-udem/fuel.git@stable
    ```
* Download the hdf5 dataset of outdoor natural images: [ourdoor_64.hdf5](http://efrosgans.eecs.berkeley.edu/iGAN/datasets/outdoor_64.zip) (1.4G), which contains 150K landscape images from MIT [Places](http://places.csail.mit.edu/) dataset. 
* Train unsupervised Blending GAN:
    ```bash
    python train_wasserstein_gan.py --data_root [Path for outdoor_64.hdf5]
    ```
* Training Curve
![](images/unsupervised_gan_result/d_loss.png)
* Samples after training

  ![](images/unsupervised_gan_result/samples.png)

## Visual results

| Mask | Copy-and-Paste | Modified-Poisson | Multi-splines | Supervised GP-GAN | Unsupervised GP-GAN |
| --- | --- | --- | --- | --- | --- |
| ![](images/result_comparison/740_mask.png) | ![](images/result_comparison/740_copy-paste.png) | ![](images/result_comparison/740_modified-poisson.png) | ![](images/result_comparison/740_multi-splines.png) | ![](images/result_comparison/740_poisson-gan-encoder.png) | ![](images/result_comparison/740_poisson-gan-wgan.png) |
| ![](images/result_comparison/2357_mask.png) | ![](images/result_comparison/2357_copy-paste.png) | ![](images/result_comparison/2357_modified-poisson.png) | ![](images/result_comparison/2357_multi-splines.png) | ![](images/result_comparison/2357_poisson-gan-encoder.png) | ![](images/result_comparison/2357_poisson-gan-wgan.png) |
| ![](images/result_comparison/1550_mask.png) | ![](images/result_comparison/1550_copy-paste.png) | ![](images/result_comparison/1550_modified-poisson.png) | ![](images/result_comparison/1550_multi-splines.png) | ![](images/result_comparison/1550_poisson-gan-encoder.png) | ![](images/result_comparison/1550_poisson-gan-wgan.png) |
| ![](images/result_comparison/1920_mask.png) | ![](images/result_comparison/1920_copy-paste.png) | ![](images/result_comparison/1920_modified-poisson.png) | ![](images/result_comparison/1920_multi-splines.png) | ![](images/result_comparison/1920_poisson-gan-encoder.png) | ![](images/result_comparison/1920_poisson-gan-wgan.png) |
| ![](images/result_comparison/1153_mask.png) | ![](images/result_comparison/1153_copy-paste.png) | ![](images/result_comparison/1153_modified-poisson.png) | ![](images/result_comparison/1153_multi-splines.png) | ![](images/result_comparison/1153_poisson-gan-encoder.png) | ![](images/result_comparison/1153_poisson-gan-wgan.png) |
