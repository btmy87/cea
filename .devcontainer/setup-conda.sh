#!/bin/bash
# Setup conda environment
# Invoked during .devcontainer postCreateCommand to we can use environment.yml

# Create conda environment
conda env create -f=environment.yml
conda init bash \
echo "conda activate cea-dev" >> ~/.bashrc

# Optionally configure intel compilers
if [ ${FC} = "ifx" ] ||  [ ${CC} = "icx" ]; then 
    conda install -y dpcpp_linux-64 ifx_linux-64; 
fi
