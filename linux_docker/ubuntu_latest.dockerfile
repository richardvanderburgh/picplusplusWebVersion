FROM ubuntu:latest

RUN apt-get update && \
apt-get install -y python3 python3-pip 

RUN pip install conan 
RUN apt-get install cmake -y

# RUN mkdir build && cd build && \
# 	conan install .. --profile=../buildUtils/linux_gcc_release --build=missing -of=. && \
# 	conan build -of . ..