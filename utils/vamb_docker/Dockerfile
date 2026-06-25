FROM pytorch/pytorch:2.7.0-cuda12.6-cudnn9-runtime
LABEL authors="Eike Matthias Wacker" \
      description="Docker image containing all requirements for recent version of gpu capable vamb for metagenome-pipeline"
      
# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PIP_NO_CACHE_DIR=false
# Install system dependencies
RUN apt-get update && apt-get install -y git
# Install vamb
RUN git clone https://github.com/RasmussenLab/vamb -b v5.0.4 && cd vamb && pip install -e . --no-cache-dir
# Make vamb utility scripts executable
ENV PATH="/workspace/vamb/src:${PATH}"
ARG PATH="/workspace/vamb/src:${PATH}"
RUN ["chmod", "+x", "/workspace/vamb/src/concatenate.py"]
RUN ["chmod", "+x", "/workspace/vamb/src/merge_aemb.py"]
RUN ["chmod", "+x", "/workspace/vamb/vamb/vambtools.py"]