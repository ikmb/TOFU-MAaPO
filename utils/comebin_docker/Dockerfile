FROM ubuntu:noble
LABEL authors="Eike Matthias Wacker" \
      description="Docker image containing all requirements for gpu capable comebin for the IKMB metagenome pipeline"

WORKDIR /app

# Copy the environment file
COPY pixi.toml /app

RUN apt-get update && apt-get install -y \
    curl \
    git \
    build-essential \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*


# Install Pixi
ENV PIXI_HOME=/opt/pixi
RUN mkdir ${PIXI_HOME} && curl -Ls https://github.com/prefix-dev/pixi/releases/download/v0.47.0/pixi-x86_64-unknown-linux-musl -o ${PIXI_HOME}/pixi
RUN chmod +x ${PIXI_HOME}/pixi
ENV PATH="${PIXI_HOME}:${PATH}"

# Install the environment
ENV CONDA_OVERRIDE_CUDA=12.6
RUN pixi install
ENV PATH="/app/.pixi/envs/default/bin:${PATH}"
# Install checkm database
RUN checkm data setRoot /app/checkm_data
RUN curl -L https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz | tar -xz -C /app/checkm_data --strip-components=1
