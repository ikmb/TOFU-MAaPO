FROM nfcore/base
LABEL authors="Eike Wacker" \
      description="Docker image containing all requirements for IKMB metagenome pipeline"

COPY environment.yml /
RUN conda install -c conda-forge mamba && conda clean -a
RUN mamba env create -f /environment.yml && mamba clean -a
#ENV PATH /opt/conda/envs/ikmb-metagenome-1.2/bin:$PATH
ENV PATH="${PATH}:/opt/conda/envs/ikmb-metagenome-1.2/bin"
