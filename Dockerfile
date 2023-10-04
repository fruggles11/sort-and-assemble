# Use an official Ubuntu as a parent image
FROM ubuntu:latest

# Set working directory
WORKDIR /scratch

# Set the maintainer label
LABEL maintainer="nrminor@wisc.edu"

# Set time zone
ENV TZ America/New_York

# Set environment variables to non-interactive (this prevents some prompts)
ENV DEBIAN_FRONTEND=non-interactive

# Run package updates and install packages
RUN apt-get update \
    && apt-get install -y \
    wget \
    default-jdk \
    unzip \
    zstd \
    build-essential \
    zlib1g-dev \
    git

# Install BBTools
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz \
    && tar -xvzf BBMap_38.90.tar.gz \
    && rm BBMap_38.90.tar.gz \
    && mv bbmap /opt/

# Update PATH for BBTools
ENV PATH="/opt/bbmap:${PATH}"

# Install seqkit
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.1.0/seqkit_linux_amd64.tar.gz \
    && tar -xvzf seqkit_linux_amd64.tar.gz \
    && mv seqkit /usr/local/bin/ \
    && rm seqkit_linux_amd64.tar.gz

# Install csvtk
RUN wget https://github.com/shenwei356/csvtk/releases/download/v0.23.0/csvtk_linux_amd64.tar.gz \
    && tar -xvzf csvtk_linux_amd64.tar.gz \
    && mv csvtk /usr/local/bin/ \
    && rm csvtk_linux_amd64.tar.gz

# Install Canu
RUN git clone https://github.com/marbl/canu.git \
    && cd canu/src \
    && make -j \
    && cd ../Linux-amd64/bin \
    && cp canu /usr/local/bin/

# Install igBLAST
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-1.17.0-x64-linux.tar.gz \
    && tar -xvzf ncbi-igblast-1.17.0-x64-linux.tar.gz \
    && mv ncbi-igblast-1.17.0 /opt/ \
    && rm ncbi-igblast-1.17.0-x64-linux.tar.gz

# Update PATH for igBLAST
ENV PATH="/opt/ncbi-igblast-1.17.0/bin:${PATH}"

# Default command to execute when starting a container from this image
CMD [ "bash" ]
