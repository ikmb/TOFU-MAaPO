FROM mambaorg/micromamba:2.0.5
LABEL authors="Eike Wacker" \
      description="Docker image containing all requirements for IKMB metagenome pipeline"

COPY environment.yml /

RUN micromamba env create -y -f /environment.yml && \
    micromamba clean -a -y

ENV PATH="/opt/conda/envs/ikmb-metagenome-1.2/bin:${PATH}"
