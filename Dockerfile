FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

RUN apt-get update && apt install -y cmake \
    libboost-system-dev libboost-thread-dev libboost-serialization-dev \
    libboost-filesystem-dev libboost-program-options-dev libboost-timer-dev \
    libsm6 \
    wget \
    openbabel \
    apt git python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python

# Install uni-dock
RUN cd opt && \
    wget https://github.com/dptech-corp/Uni-Dock/releases/download/1.0.0/unidock && \
    chmod +x unidock

# Download AFDR tooling
RUN cd opt && \
    wget https://ccsb.scripps.edu/adfr/download/1038/ -O adfr.tar.gz && \
    tar -xzvf adfr.tar.gz && \
    rm adfr.tar.gz && \
    cd ADFRsuite_x86_64Linux_1.0 && \
    echo "Y" | ./install.sh -d afdr -c 0

# Install python dependencies
RUN pip install duckdb pandas hydra-core tqdm rdkit --upgrade

# Ensure binaries are in path
ENV PATH="/opt:${PATH}"
ENV PATH="/opt/ADFRsuite_x86_64Linux_1.0/afdr/bin:${PATH}"

# Copy source folder into app folder
COPY ./src /app/src

# Create an output directory
RUN mkdir ./data

# Move into app folder
WORKDIR /app

# Run the app
CMD ["python3", "-m", "src.unidock.run"]