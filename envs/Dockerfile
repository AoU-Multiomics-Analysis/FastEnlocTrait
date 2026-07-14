FROM mambaorg/micromamba:1.5.3

ENV MAMBA_DOCKERFILE_ACTIVATE=1

ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    conda-forge::r-tidyverse \
    conda-forge::r-data.table \
    conda-forge::r-optparse \
    conda-forge::r-ggh4x \
    conda-forge::r-arrow \
    conda-forge::r-r.utils \
    conda-forge::r-bedr \
    conda-forge::r-readxl \
    bioconda::htslib \
    conda-forge::google-cloud-sdk \
    conda-forge::curl \
    conda-forge::gzip \
    conda-forge::git \
    conda-forge::make \
    conda-forge::gcc \
    conda-forge::gxx \
    conda-forge::gsl \
    conda-forge::zlib \
    && micromamba clean -a -y

USER root

RUN git clone https://github.com/xqwen/fastenloc.git /opt/fastenloc && \
    cd /opt/fastenloc/src && \
    make && \
    cp fastenloc /usr/local/bin/fastenloc && \
    chmod +x /usr/local/bin/fastenloc

USER mambauser

WORKDIR /home/mambauser

COPY scripts/* .

CMD ["bash"]
