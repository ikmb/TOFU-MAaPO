FROM nfcore/base
LABEL authors="Eike Wacker" \
      description="Docker image containing all requirements for IKMB metagenome pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ikmb-metagenome-1.2/bin::$PATH
