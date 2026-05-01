FROM ubuntu:24.04

# Avoid interactive prompts during package install
ENV DEBIAN_FRONTEND=noninteractive

# Install build tools and Eigen3
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libeigen3-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy project into container
WORKDIR /project
COPY . .

# Build
RUN mkdir -p build && cd build && cmake .. && cmake --build .

# Default command: show usage
CMD ["./build/vibrational_frequency"]
