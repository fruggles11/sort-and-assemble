# Use an official Ubuntu as a parent image
FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set the maintainer label
LABEL maintainer="nrminor@wisc.edu"

# Set time zone
ENV TZ America/New_York

# Set environment variables to non-interactive (this prevents some prompts)
ENV DEBIAN_FRONTEND=non-interactive

# Run package updates and install packages
# Use `dpkg -l > apt-get.lock` inside a container to print package versions
COPY resources/config/apt-get.lock /apt-get-tmp/package_list.txt
RUN apt-get update && \
    awk '/^ii/ { printf("apt-get install -y %s=%s\n", $2, $3) }' /apt-get-tmp/package_list.txt \
    > /apt-get-tmp/install_packages.sh && \
    chmod +x /apt-get-tmp/install_packages.sh && \
    /apt-get-tmp/install_packages.sh && \
    rm -rf /apt-get-tmp && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# build go package for downloading files
COPY src/ /src/
RUN cd /src && \
    go build goDownloadFiles.go && \
    chmod +x goDownloadFiles
ENV PATH="${PATH}:/src"

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
RUN curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.tar.xz --output /opt/canu-2.2.tar.xz && \
    cd /opt && \
    tar -xJf canu-2.2.tar.xz && \
    rm canu-2.2.tar.xz && \
    cd canu-2.2/src && \
    make -j 8
ENV PATH=$PATH:/opt/canu-2.2/build/bin

# Install igBLAST
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.9.0/ncbi-igblast-1.9.0-x64-linux.tar.gz \
    && tar -xvzf ncbi-igblast-1.9.0-x64-linux.tar.gz \
    && mv ncbi-igblast-1.9.0 /opt/ \
    && rm ncbi-igblast-1.9.0-x64-linux.tar.gz

# Update PATH for igBLAST
ENV PATH="/opt/ncbi-igblast-1.9.0/bin:${PATH}"

# Default command to execute when starting a container from this image
CMD [ "bash" ]
