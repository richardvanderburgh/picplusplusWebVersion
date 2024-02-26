FROM ubuntu:latest

RUN apt-get update && \
apt-get install -y python3 python3-pip 

RUN pip install conan 
RUN apt-get install cmake -y
