# Dockerfile
FROM mambaorg/micromamba:1.5.8

# 1. Copy the conda environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.yaml.txt /tmp/env.yaml

# 2. Install dependencies (fastp, bowtie2, python, etc.)
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

# 3. Copy your specific source code into the container
#    (Assuming you moved your code to src/ structure, otherwise adjust path)
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# 4. Install your python tool inside the container
WORKDIR /app
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install .

# 5. Set the entrypoint so running the container runs your tool
ENTRYPOINT ["trnagraph"]