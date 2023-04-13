FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Install general tools
RUN apt-get update && apt-get install -y \
    build-essential \
    libbz2-dev \
    wget \
    zlib1g-dev

# Install TrioDenovo
RUN mkdir triodenovo &&\
    cd triodenovo && \
    wget https://genome.sph.umich.edu/w/images/f/f4/Triodenovo.0.06.tar.gz &&\
    tar -xf Triodenovo.0.06.tar.gz

COPY Argument.h /triodenovo/triodenovo.0.06/base/
    
RUN cd triodenovo/triodenovo.0.06/ &&\
    make &&\
    rm /triodenovo/Triodenovo.0.06.tar.gz

ENV TRIODENOVO="/triodenovo/triodenovo.0.06/bin/triodenovo"