FROM continuumio/miniconda3

COPY conda_linux64.yaml /tmp/conda_linux64.yaml
RUN conda update -n base -c defaults conda
RUN conda env create -f tmp/conda_linux64.yaml 
RUN echo "source activate mgefinder" >> ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
