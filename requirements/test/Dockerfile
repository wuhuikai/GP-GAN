FROM chainer/chainer
MAINTAINER romusters@gmail.com

# Python3.5
RUN apt-get update -y
RUN apt-get install -y software-properties-common apt-utils
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update -y
RUN apt-get install -y python3.5 python3.5-dev python3-pip git

RUN git clone https://github.com/wuhuikai/GP-GAN.git
WORKDIR GP-GAN
RUN pip3 install --upgrade pip
RUN pip3 install -r requirements/test/requirements.txt
