# Chloe (https://github.com/ian-small/Chloe.jl) for Nextflow/Singularity.
# The upstream image keeps its Julia packages in /root/.julia with /root at mode 700,
# so a non-root user (Singularity, docker -u) cannot load them, and it lacks the
# `ps` command Nextflow needs for task metrics. This layer fixes both.
FROM docker.io/willmundlab/chloe:0.1.17

RUN chmod 755 /root \
    && apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

# First depot entry is writable scratch for (re)compilation; the bundled one follows.
ENV JULIA_DEPOT_PATH=/tmp/julia:/root/.julia
