# Use Ubuntu as the base image
FROM ubuntu:latest

# Install necessary packages and dependencies
RUN apt-get update && apt-get install -y \
    wget \
    software-properties-common \
    curl \
    git \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /SpecImmune

# Define SHELL
SHELL ["/bin/bash", "-c"] 

# Copy code and scripts into the container
COPY . /SpecImmune

# 1. Download and install Miniforge
RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh -O /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh

# Add Conda to the PATH
ENV PATH="/opt/conda/bin:$PATH"

# Create the env from conda-forge/bioconda
RUN conda config --system --set channel_priority strict \
    && mamba env create -p /opt/conda/envs/SpecImmune -f /SpecImmune/docker/environment.docker.yml \
    && /opt/conda/envs/SpecImmune/bin/pip install --no-deps dysgu==1.6.2 \
    && conda clean --all -y

# "Activate" the env for every later step and for the running container
ENV PATH=/opt/conda/envs/SpecImmune/bin:$PATH
ENV CONDA_DEFAULT_ENV=SpecImmune