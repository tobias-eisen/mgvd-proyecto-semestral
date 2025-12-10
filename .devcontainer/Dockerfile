FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# System setup and dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential make gcc g++ libc6-dev zlib1g-dev ca-certificates \
    git curl wget vim less python3 python3-pip pkg-config file bzip2 locales \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /mgvd

# perl configuration
RUN locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 && \
    rm -rf /var/lib/apt/lists/*

ENV LANG=en_US.UTF-8 \
    LANGUAGE=en_US:en \
    LC_ALL=en_US.UTF-8

# Install Miniforge (ARM variant) and setup channels
ENV CONDA_DIR=/opt/conda
RUN curl -L -o /tmp/miniforge.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh && \
    bash /tmp/miniforge.sh -b -p $CONDA_DIR && \
    rm /tmp/miniforge.sh && \
    $CONDA_DIR/bin/conda clean -afy

# Add conda to PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# Setup bioconda + conda-forge channels and install KrakenUniq
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y krakenuniq && \
    conda install numpy scikit-learn matplotlib && \
    conda clean -afy

RUN conda init bash

CMD ["/bin/bash"]
