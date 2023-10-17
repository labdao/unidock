FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

RUN apt-get update && apt install -y cmake \
    libboost-system-dev libboost-thread-dev libboost-serialization-dev \
    libboost-filesystem-dev libboost-program-options-dev libboost-timer-dev \
    python3 git python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python

# Install python dependencies
RUN pip install click

# Install uni-dock
RUN cd opt && \
    git clone https://github.com/dptech-corp/Uni-Dock.git && \
    cd /opt/Uni-Dock/unidock && \
    cmake -B build && \
    cmake --build build -j`nprocs` && \
    cmake --install build && \
    rm -r /opt/Uni-Dock

# Copy source files into app folder
COPY . /app

# Move into app folder
WORKDIR /app