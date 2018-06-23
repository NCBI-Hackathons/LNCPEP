FROM ubuntu:16.04

MAINTAINER Stephen M. Kelly

# ~~~~~ BASIC SETUP ~~~~~ #
RUN apt-get update && \
apt-get install -y wget \
bzip2


# ~~~~~ MINICONDA ~~~~~ #
RUN wget https://repo.continuum.io/miniconda/Miniconda3-4.4.10-Linux-x86_64.sh && \
bash Miniconda3-4.4.10-Linux-x86_64.sh -b -p /conda && \
rm -f Miniconda3-4.4.10-Linux-x86_64.sh
ENV PATH="/conda/bin:${PATH}"
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda


# ~~~~~ TOOLS ~~~~~ #
RUN conda install -y -c bioconda \
sra-tools==2.9.0 \
hisat2==2.0.5 \
samtools==1.3.1 \
stringtie==1.2.4 \
sambamba==0.6.6 \
cufflinks==2.2.1

