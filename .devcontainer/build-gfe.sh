#!/bin/bash
# Build GFE from source
# This is invoked during .devcontainer postCreateCommand
# Because we needed the functioning conda environment
# which is also configured in the post-create command

# Checkout GFE
git clone https://github.com/Goddard-Fortran-Ecosystem/GFE extern/gfe
git -C extern/gfe submodule update --init 

# Configure and build
cmake -S extern/gfe -B ${BUILDDIR}/extern/gfe -G Ninja\
        -DCMAKE_C_COMPILER=${CC} \
        -DCMAKE_CXX_COMPILER=${CXX} \
        -DCMAKE_Fortran_COMPILER=${FC} \
        -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} \
        -DCMAKE_BUILD_TYPE=Release \
        -DSKIP_MPI=YES \
        -DENABLE_MPI=OFF \
        -DSKIP_OPENMP=YES \
        -DSKIP_FHAMCREST=YES \
        -DSKIP_ESMF=YES \
        -DSKIP_ROBUST=YES
        
cmake --build ${BUILDDIR}/extern/gfe --target install