FROM opencfd/openfoam-dev:2412

# Install any additional dependencies if needed
USER root
RUN apt-get update && apt-get install -y \
    build-essential \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# Create a non-root user
RUN useradd -m -s /bin/bash -G sudo ofuser && \
    echo "ofuser ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

# Switch to non-root user
USER ofuser
WORKDIR /home/ofuser

# Copy the clotFoam source code
COPY --chown=ofuser:ofuser clotFoam /home/ofuser/clotFoam
COPY --chown=ofuser:ofuser tutorials /home/ofuser/tutorials
COPY --chown=ofuser:ofuser docker/run-test.sh /home/ofuser/docker/run-test.sh
COPY --chown=ofuser:ofuser docker/test.sh /home/ofuser/docker/test.sh

# Make scripts executable
RUN chmod +x /home/ofuser/docker/run-test.sh /home/ofuser/docker/test.sh

# Create OpenFOAM user directory structure
RUN mkdir -p /home/ofuser/OpenFOAM/ofuser-v2412/platforms/linuxARM64GccDPInt32Opt/bin && \
    mkdir -p /home/ofuser/OpenFOAM/ofuser-v2412/platforms/linuxARM64GccDPInt32Opt/lib

# Set environment variables
ENV WM_PROJECT_USER_DIR=/home/ofuser/OpenFOAM/ofuser-v2412
ENV FOAM_USER_APPBIN=/home/ofuser/OpenFOAM/ofuser-v2412/platforms/linuxARM64GccDPInt32Opt/bin
ENV FOAM_USER_LIBBIN=/home/ofuser/OpenFOAM/ofuser-v2412/platforms/linuxARM64GccDPInt32Opt/lib
