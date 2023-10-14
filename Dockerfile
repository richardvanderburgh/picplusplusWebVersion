# Use an official Python runtime as a parent image
# FROM python:3.9-slim
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

EXPOSE 8000

# Install Python (Python 3 in this case) and required dependencies
RUN apt-get update && \
apt-get install -y python3 python3-pip vim
# RUN apt-get install curl

# Install Conan using pip
RUN pip3 install conan==1.54 \
&& pip3 install Django \
&& pip3 install django-cors-headers \
&& pip3 install gunicorn==20.1.0 \
&& pip install whitenoise


# Install CMake
RUN apt-get install -y cmake

# Set a working directory (optional)
WORKDIR /app
# # Copy the current directory contents into the container at /app
COPY . /app

RUN mkdir build \
&& cd build \ 
&& conan install .. --profile=../linux_profile --build=missing \
&& cmake -B . -S .. \ 
&& make \
&& mv /app/build/bin/PIC++Main /

# RUN mv /app/build/bin/PIC++Main

# Expose port 8000
EXPOSE 8000

# Start Gunicorn to serve your Django app
CMD ["/app/start.sh"]
#CMD ["gunicorn", "djangodocker.wsgi:application", "--bind", "0.0.0.0:$PORT"]
