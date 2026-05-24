#!/bin/bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc

# Override user directories
export WM_PROJECT_USER_DIR=/home/ofuser/OpenFOAM/ofuser-v2412
export FOAM_USER_APPBIN=/home/ofuser/OpenFOAM/ofuser-v2412/platforms/linuxARM64GccDPInt32Opt/bin
export FOAM_USER_LIBBIN=/home/ofuser/OpenFOAM/ofuser-v2412/platforms/linuxARM64GccDPInt32Opt/lib
export PATH=$FOAM_USER_APPBIN:$PATH

exec /home/ofuser/docker/test.sh
