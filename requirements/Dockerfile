FROM chainer/chainer
MAINTAINER romusters@gmail.com

RUN apt-get update

# Python3.5
RUN apt-get -y install software-properties-common apt-utils nano
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt-get update
RUN apt-get install -y python3.5 python3.5-dev

RUN apt-get -y install python3-pip git

COPY . .
RUN pip3 install -r requirements/requirements.txt
