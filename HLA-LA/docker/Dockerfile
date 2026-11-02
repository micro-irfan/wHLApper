# HLA*LA image for wHLApper (processes with label 'hlala')
#
# Contains: HLA*LA 1.0.3, samtools, bwa, picard, minimap2 (for MINIMAP2_ALIGN)
# The HLA*LA graph (PRG_MHC_GRCh38_withIMGT) is NOT included - pass it with --graph_dir.
#
# Build:
#   cd HLA-LA && docker build -t hla-la:1.0.3 docker
#
# Adapted from reference/Dockerfile, without the private minimap2 base image.

FROM ubuntu:22.04

ARG SAMTOOLS_VERSION=1.21
ARG MINIMAP2_VERSION=2.28
ARG PICARD_VERSION=2.27.5
ARG BOOST_VERSION=1.80.0
ARG BAMTOOLS_VERSION=2.5.2
# Commit tied to HLA*LA v1.0.3
ARG HLALA_COMMIT=fb8ae1cde7783d6cf441a80319e943f63f3b6276

ENV DEBIAN_FRONTEND=noninteractive

# Build tools, libraries, Java for picard, bwa, and HLA*LA's Perl dependencies.
# procps provides `ps`, which Nextflow needs inside the container.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        git \
        wget \
        bzip2 \
        procps \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libicu-dev \
        libjsoncpp25 \
        openjdk-17-jre-headless \
        bwa \
        perl \
        cpanminus \
        bioperl \
        libbio-db-hts-perl \
        libbio-featureio-perl \
        libbio-samtools-perl \
        libdata-dump-perl \
        liblist-moreutils-perl \
        libmath-basecalc-perl \
        libparallel-forkmanager-perl \
        libtext-levenshtein-perl \
    && cpanm --notest Array::Utils \
    && perl -MArray::Utils -e ';' \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /root/.cpanm

WORKDIR /opt2

# samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && cd samtools-${SAMTOOLS_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j"$(nproc)" \
    && make install \
    && cd .. \
    && rm -rf samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.tar.bz2

# minimap2 (prebuilt x86_64 binary)
RUN wget -q https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 \
    && tar xjf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 \
    && cp minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2 /usr/local/bin/ \
    && rm -rf minimap2-${MINIMAP2_VERSION}_x64-linux*

# picard (2.x keeps the legacy I=/F= syntax HLA*LA uses for SamToFastq)
RUN mkdir -p /opt2/picard \
    && wget -q https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar -O /opt2/picard/picard.jar
ENV PICARD=/opt2/picard/picard.jar

# Boost, built from source: HLA*LA's makefile needs <prefix>/include and <prefix>/lib
RUN BOOST_US=$(echo ${BOOST_VERSION} | tr . _) \
    && wget -q https://archives.boost.io/release/${BOOST_VERSION}/source/boost_${BOOST_US}.tar.bz2 \
    && tar xjf boost_${BOOST_US}.tar.bz2 \
    && cd boost_${BOOST_US} \
    && ./bootstrap.sh --prefix=/opt2/libboost \
    && ./b2 -q -j"$(nproc)" install --prefix=/opt2/libboost \
    && cd .. \
    && rm -rf boost_${BOOST_US} boost_${BOOST_US}.tar.bz2

# bamtools (shared libraries)
RUN wget -q https://github.com/pezmaster31/bamtools/archive/refs/tags/v${BAMTOOLS_VERSION}.tar.gz \
    && tar xzf v${BAMTOOLS_VERSION}.tar.gz \
    && mkdir bamtools-${BAMTOOLS_VERSION}/build \
    && cd bamtools-${BAMTOOLS_VERSION}/build \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt2/libbamtools \
        -DCMAKE_INSTALL_INCLUDEDIR=/opt2/libbamtools/include \
        -DCMAKE_INSTALL_LIBDIR=/opt2/libbamtools/lib \
        -DBUILD_SHARED_LIBS=TRUE \
        .. \
    && make -j"$(nproc)" install \
    && cd ../.. \
    && rm -rf bamtools-${BAMTOOLS_VERSION} v${BAMTOOLS_VERSION}.tar.gz \
    && ln -sf /opt2/libbamtools/lib /opt2/libbamtools/lib64

ENV PATH="/opt2/libboost/bin:/opt2/libbamtools/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt2/libboost/lib:/opt2/libbamtools/lib"

# HLA*LA: compile, point paths.ini at the tools above, and test the binary.
# workingDir is left empty because wHLApper always passes --workingDir.
RUN mkdir -p /opt2/hla-la/HLA-LA/bin \
        /opt2/hla-la/HLA-LA/src \
        /opt2/hla-la/HLA-LA/obj \
        /opt2/hla-la/HLA-LA/temp \
        /opt2/hla-la/HLA-LA/working \
        /opt2/hla-la/HLA-LA/graphs \
    && cd /opt2/hla-la/HLA-LA/src \
    && git clone https://github.com/DiltheyLab/HLA-LA.git . \
    && git checkout ${HLALA_COMMIT} \
    && make -j"$(nproc)" all \
        BOOST_PATH=/opt2/libboost \
        BAMTOOLS_PATH=/opt2/libbamtools \
    && sed -i "s:^\(picard_sam2fastq_bin=\).*:\1${PICARD}:" paths.ini \
    && sed -i "s:^\(samtools_bin=\).*:\1/usr/local/bin/samtools:" paths.ini \
    && sed -i "s:^\(bwa_bin=\).*:\1/usr/bin/bwa:" paths.ini \
    && sed -i "s:^\(workingDir=\).*:workingDir=:" paths.ini \
    && sed -i "s:^\(workingDir_HLA_ASM=\).*:workingDir_HLA_ASM=:" paths.ini \
    && chmod +x *.pl \
    && rm -rf .git /opt2/hla-la/HLA-LA/obj/* \
    && /opt2/hla-la/HLA-LA/bin/HLA-LA --action testBinary

ENV PATH="/opt2/hla-la/HLA-LA/bin:/opt2/hla-la/HLA-LA/src:${PATH}"

# Quick check that every tool the pipeline calls is on PATH
RUN samtools --version | head -n1 \
    && minimap2 --version \
    && (bwa 2>&1 | grep Version) \
    && java -jar ${PICARD} SamToFastq --version 2>&1 | tail -n1 \
    && which HLA-LA.pl

WORKDIR /data
