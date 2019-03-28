FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
  locales \
  sudo \
  git \
  vim \
  python3.5 \
  python3-pip \
  python3-dev

RUN pip3 install --upgrade pip
RUN python3 -V
COPY ./requirements.txt /requirements.txt
RUN pip3 install -r /requirements.txt

ENV PYTHONIOENCODING=utf-8
RUN locale-gen en_US.UTF-8
ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'
RUN echo "ja_JP UTF-8" > /etc/locale.gen


WORKDIR /work

ENTRYPOINT ["/bin/bash"]
