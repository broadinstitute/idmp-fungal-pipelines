FROM ubuntu:22.04

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    snp-sites \
    samtools \
    tabix \
    curl \
    && apt-get clean
