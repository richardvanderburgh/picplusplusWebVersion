FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
apt-get install -y python3 python3-pip vim

RUN pip3 install conan==1.54

RUN apt-get install -y cmake

WORKDIR /app
COPY . /app

RUN mkdir buildDocker \
&& cd buildDocker \ 
&& conan install .. --profile=../buildUtils/linux_profile --build=missing \
&& cmake -B . -S .. \ 
&& make \
&& mv /app/build/bin/PIC++Main /