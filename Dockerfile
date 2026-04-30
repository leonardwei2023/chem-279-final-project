FROM ubuntu:22.04

# Avoiding prompts:

ENV DEBIAN_FRONTEND=noninteractive

# Installing dependencies:

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

# Setting working directory:

WORKDIR /app

# Copying project files:

COPY . .

# Building:

RUN mkdir build && cd build && cmake .. && make

# Default command:

CMD ["bash"]
