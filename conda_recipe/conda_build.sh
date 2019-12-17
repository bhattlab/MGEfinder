#!/bin/bash
conda-build mgefinder/ -c bioconda -c conda-forge
mkdir -p outputdir/linux-64
cp /home/mdurrant/miniconda3/envs/mgefinder/conda-bld/linux-64/mgefinder-1.0.0-py36_0.tar.bz2 outputdir/linux-64
conda convert --platform all /home/mdurrant/miniconda3/envs/mgefinder/conda-bld/linux-64/mgefinder-1.0.0-py36_0.tar.bz2 -o outputdir/
