FROM nvidia/cuda:10.0-cudnn7-runtime-ubuntu18.04
MAINTAINER r.musters@anchormen.nl

RUN apt-get update

# Python3.5
RUN apt-get -y install software-properties-common apt-utils
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update
RUN apt-get install -y python3.5 python3.5-dev

RUN apt-get -y install python3-pip git

#RUN python3.5 -m pip uninstall pip && apt install python3-pip --reinstall
#ADD requirements.txt .

RUN git clone https://github.com/wuhuikai/GP-GAN.git
RUN python3.5 -m pip install -r requirements.txt


