FROM ubuntu:latest

LABEL maintainer="Andre Fonseca" \
    description="scAligners - Container with Cellranger"

ARG CELLRANGER_VERSION='7.1.0'

## Please change the link considering the 10x Genomics End User Software License Agreement
ARG CELLRANGER_URL='https://cf.10xgenomics.com/releases/cell-exp/cellranger-${CELLRANGER_VERSION}.tar.gz'
ENV PATH /opt/cellranger-${CELLRANGER_VERSION}:$PATH

# Install basic requirments
RUN apt-get update && \
    apt-get install -y wget

# CellRanger binaries
RUN wget -O cellranger-${CELLRANGER_VERSION}.tar.gz ${CELLRANGER_URL} \
    && tar xzf cellranger-${CELLRANGER_VERSION}.tar.gz \
    && rm -rf cellranger-${CELLRANGER_VERSION}.tar.gz \
    && mv cellranger-${CELLRANGER_VERSION} /opt/

WORKDIR /opt

CMD  ["cellranger" ]