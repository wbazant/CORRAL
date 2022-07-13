FROM ubuntu:22.04

LABEL maintainer="rdemko2332@gmail.com"


#Installing Software
RUN apt-get update && \
    apt-get install -y \
    wget \
    bowtie2=2.4.4-1 \
    samtools=1.13-4 \
    bzip2 \
    sra-toolkit \
    python3-pip \
    cpanminus \
  && rm -rf /var/lib/apt/lists/*
  
RUN pip install marker_alignments && cpanm List::MoreUtils

COPY /bin/* /usr/bin/

RUN cd /usr/bin/ && chmod +x *

WORKDIR /work