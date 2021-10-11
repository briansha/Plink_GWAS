# Filename: Dockerfile
# Author: brian.sharber@vumc.org
# Used for:
# - briansha/plink2:polygenic
# - briansha/plink2:terra

FROM ubuntu:18.04

WORKDIR /cromwell_root/

RUN apt-get update && apt-get install -y --no-install-recommends \
      curl \
      zip \
      git \
      unzip \
      gzip \
      g++ \
      make \
      python3 \
      python3-pip \
      gfortran \
      zlib1g-dev \
      libgfortran4 \
      liblapacke-dev \
      libopenblas-dev \
      libbz2-dev \
      liblzma-dev

# Python packages
RUN pip3 install scipy
RUN python3 -m pip install h5py
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3-pandas

# PLINK2
ADD http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip .
RUN unzip plink2_linux_avx2.zip && rm plink2_linux_avx2.zip
RUN cp plink2 /usr/local/bin

# HTSLib
ADD https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 .
RUN tar -xvf htslib-1.12.tar.bz2 && rm htslib-1.12.tar.bz2
WORKDIR /cromwell_root/htslib-1.12
RUN ./configure --prefix=/usr/local
RUN make
RUN make install
WORKDIR /cromwell_root/
RUN rm -r htslib-1.12
