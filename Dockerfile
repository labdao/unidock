FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

RUN apt-get update && apt install -y cmake \
    libboost-system-dev libboost-thread-dev libboost-serialization-dev \
    libboost-filesystem-dev libboost-program-options-dev libboost-timer-dev \
    libsm6 \
    wget \
    openbabel \
    python3 git python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python

# Install uni-dock
RUN git clone https://github.com/dptech-corp/Uni-Dock.git && \
    cd Uni-Dock/unidock && \
    cmake -B build && \
    cmake --build build -j4 && \
    cmake --install build && \
    cd ../unidock_tools && \
    python setup.py install

# Install python dependencies
RUN pip install duckdb pandas hydra-core tqdm --upgrade

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