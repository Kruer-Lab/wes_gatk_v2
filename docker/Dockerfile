#FROM ubuntu:22.04
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Install general tools
RUN apt-get update && apt-get install -y \
    build-essential \
    bzip2 \
    less \
    libbz2-dev \
    liblzma-dev \
    default-jdk \
    pkg-config \
    python-is-python3 \
    unzip \
    wget \
    zlib1g-dev

# VCFtools v0.1.16
RUN mkdir vcftools &&\
    cd vcftools &&\
    wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz &&\
    tar -xf vcftools-0.1.16.tar.gz &&\
    rm vcftools-0.1.16.tar.gz &&\
    cd vcftools-0.1.16 &&\
    ./configure &&\
    make &&\
    make install

# BWA version 0.7.17 (r1188)
RUN mkdir bwa &&\
    cd bwa &&\
    wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 &&\
    tar -xf bwa-0.7.17.tar.bz2 &&\
    rm bwa-0.7.17.tar.bz2 &&\
    cd bwa-0.7.17 &&\
    make CC='gcc -fcommon'

ENV PATH="$PATH:/bwa/bwa-0.7.17"

# Samtools version 1.16.1
RUN mkdir samtools &&\
    cd samtools &&\
    wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 &&\
    tar -xf samtools-1.16.1.tar.bz2 &&\
    rm samtools-1.16.1.tar.bz2 &&\
    cd samtools-1.16.1 &&\
    ./configure --without-curses --prefix=/usr/local &&\
    make &&\
    make install

# Picard version 2.27.5
RUN mkdir picard &&\
    cd picard &&\
    wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar

ENV PICARD="/picard/picard.jar"

# GATK version 4.2.6.1
RUN mkdir gatk &&\
    cd gatk &&\
    wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip &&\
    unzip gatk-4.2.6.1.zip &&\
    rm gatk-4.2.6.1.zip

ENV PATH="$PATH:/gatk/gatk-4.2.6.1"

# IGV version 4.2.6.1
RUN mkdir igvtools &&\
    cd igvtools &&\
    wget https://data.broadinstitute.org/igv/projects/downloads/2.14/IGV_2.14.1.zip &&\
    unzip IGV_2.14.1.zip &&\
    rm IGV_2.14.1.zip

ENV PATH="$PATH:/igvtools/IGV_2.14.1"

# Bedtools version 2.30.0
RUN mkdir bedtools &&\
    cd bedtools &&\
    wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary &&\
    mv bedtools.static.binary bedtools &&\
    chmod a+x bedtools

ENV PATH="$PATH:/bedtools"

# Install Annovar (Date:2020-06-07)
RUN mkdir annovar &&\
    cd annovar &&\
    wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz &&\
    tar -xf annovar.latest.tar.gz &&\
    rm annovar.latest.tar.gz

ENV ANNOVAR="/annovar/annovar/"

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