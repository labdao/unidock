FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

RUN apt-get update && apt install -y cmake \
    libboost-system-dev libboost-thread-dev libboost-serialization-dev \
    libboost-filesystem-dev libboost-program-options-dev libboost-timer-dev \
    libsm6 \
    wget \
    openbabel \
    python3 git python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python

# Install python dependencies
RUN pip install click duckdb

# Install uni-dock
RUN cd opt && \
    git clone https://github.com/dptech-corp/Uni-Dock.git && \
    cd /opt/Uni-Dock/unidock && \
    cmake -B build && \
    cmake --build build -j`nprocs` && \
    cmake --install build && \
    rm -r /opt/Uni-Dock

# Download AFDR tooling
RUN cd opt && \
    wget https://ccsb.scripps.edu/adfr/download/1038/ -O adfr.tar.gz && \
    tar -xzvf adfr.tar.gz && \
    rm adfr.tar.gz && \
    cd ADFRsuite_x86_64Linux_1.0 && \
    echo "Y" | ./install.sh -d afdr -c 0

ENV PATH="/opt/ADFRsuite_x86_64Linux_1.0/afdr/bin:${PATH}"

# Copy source files into app folder
COPY ./src /app

# Move into app folder
WORKDIR /app
