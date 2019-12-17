Bootstrap: docker

From: continuumio/miniconda3

%environment
    PATH=/opt/conda/envs/mgefinder/bin:$PATH
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8

%post
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
    echo "source activate mgefinder" > ~/.bashrc
    git clone https://github.com/bhattlab/MGEfinder.git
    /opt/conda/bin/conda env create -f MGEfinder/env/conda_linux64.yaml
    rm -rf MGEfinder

%runscript
    exec "$@"
