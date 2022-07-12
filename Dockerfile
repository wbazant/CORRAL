FROM ubuntu:22.04

MAINTAINER rdemko2332@gmail.com

RUN apt-get -qq update --fix-missing

#Installing Software
RUN apt-get install -y \
  wget \
  bowtie2=2.4.4-1 \
  samtools=1.13-4 \
  bzip2 \
  sra-toolkit \
  python3-pip \
  cpanminus
  
RUN pip install marker_alignments && cpanm List::MoreUtils