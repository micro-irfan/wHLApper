############# OVERVIEW #############
# Use an official Ubuntu as a parent image
FROM 682905754632.dkr.ecr.ap-southeast-1.amazonaws.com/nala/minimap2:2.22-r1101

# Set environment variables for non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive

# Update and install necessary packages
RUN apt-get update && \
    apt-get install -y build-essential wget cmake protobuf-compiler git zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev

# Base image already contains:
# g++ (support for C++11), libboost (apt-get version not compatible with python3.10), 
# zlib, samtools (1.17), picard
# @INSTALLS:
#  • bamtools/2.5.2  ✓
#  • bwa/0.7.17-6    ✓
#  • libboost/1.80.0 ✓
#  • HLA-LA/1.0.3    ✓
#  • bamtools/2.5.2  ✓
WORKDIR /opt2

# Install Samtools 1.17
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar xjf samtools-1.17.tar.bz2 && \
    rm samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.17

# Create the directory for Picard JAR
# RUN mkdir -p /usr/src/app/

# Install Picard
# RUN wget https://github.com/broadinstitute/picard/releases/download/2.26.1/picard.jar -O /usr/src/app/picard.jar

# Install picard/2.27.5,
# requires java8 or 1.8
RUN mkdir -p /opt2/picard/2.27.5/ \
    && wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar -P /opt2/picard/2.27.5/
ENV PICARDJARPATH=/opt2/picard/2.27.5/ \
    PICARD_JARPATH=/opt2/picard/2.27.5/ \
    PICARDJAR=/opt2/picard/2.27.5/picard.jar \
    PICARD_JAR=/opt2/picard/2.27.5/picard.jar \
    PICARD=/opt2/picard/2.27.5/picard.jar

# Install bwa/0.7.17-6 for HLA-LA,
# perl modules for HLA-LA,
# libicu-dev for boost
RUN apt-get update \
    && apt-get -y upgrade \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
        autotools-dev \
        bwa \
        bioperl \
        libicu-dev \
        libbio-db-hts-perl \
        libbio-featureio-perl \
        libbio-samtools-perl \
        libdata-dump-perl \
        libjsoncpp25 \
        liblist-moreutils-perl \
        libmath-basecalc-perl \
        libparallel-forkmanager-perl \
        libtext-levenshtein-perl \
    && cpan Array::Utils \
    && perl -MArray::Utils -e ';' \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


############## INSTALL #############
# Install C++ Boost/1.80.0 library,
# required for HLA-LA, apt install
# of libboost libraries are missing
# include/ and install/ directories 
# needed to make the HLA-LA exe, 
# also, versions of boost<1.80.0 
# are NOT compatiable with python3.10,
# for more information see: 
# https://www.boost.org/doc/libs/1_80_0/more/getting_started/unix-variants.html
# http://www.boost.org/build/doc/html/index.html
# https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=1000506
# https://github.com/boostorg/python/pull/344
RUN mkdir -p /opt2/lib/ \
    && wget "http://downloads.sourceforge.net/project/boost/boost/1.80.0/boost_1_80_0.tar.bz2" -P /opt2/lib/ \
    && tar -xjf /opt2/lib/boost_1_80_0.tar.bz2 -C /opt2/lib/ \
    && rm /opt2/lib/boost_1_80_0.tar.bz2 \
    && cd /opt2/lib/boost_1_80_0 \
    && ./bootstrap.sh \
    && ./b2 -q install --prefix=/opt2/libboost 
ENV PATH="${PATH}:/opt2/libboost/bin" \
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt2/libboost/lib"

# Install bamtools/2.5.2,
# for more info see:
# https://github.com/pezmaster31/bamtools
# https://github.com/pezmaster31/bamtools/wiki/Building-and-installing
RUN mkdir -p /opt2/bamtools/ \
    && wget https://github.com/pezmaster31/bamtools/archive/refs/tags/v2.5.2.tar.gz -P /opt2/bamtools/ \
    && tar -vxzf /opt2/bamtools/v2.5.2.tar.gz -C /opt2/bamtools/ \
    && rm /opt2/bamtools/v2.5.2.tar.gz \
    && mkdir -p /opt2/bamtools/bamtools-2.5.2/build \
    && cd /opt2/bamtools/bamtools-2.5.2/build \
    && cmake \
        -DCMAKE_INSTALL_PREFIX=/opt2/libbamtools \
        -DCMAKE_INSTALL_INCLUDEDIR=/opt2/libbamtools/include \
        -DCMAKE_INSTALL_LIBDIR=/opt2/libbamtools/lib \
        -DBUILD_SHARED_LIBS=TRUE \
        .. \
    && make install \
    && rm -rf /opt2/bamtools/bamtools-2.5.2/build
ENV PATH="${PATH}:/opt2/libbamtools/bin" \
    LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/opt2/libbamtools/lib"

# Install HLA-LA/1.0.3,
# for more info see:
# https://github.com/DiltheyLab/HLA-LA
# https://github.com/DiltheyLab/HLA-LA#installing-hlala
# Checkout commit tied to v1.0.3
# and test that the binary works!
RUN mkdir -p /opt2/hla-la/1.0.3/ \
        /opt2/hla-la/1.0.3/HLA-LA \
        /opt2/hla-la/1.0.3/HLA-LA/bin \
        /opt2/hla-la/1.0.3/HLA-LA/src \
        /opt2/hla-la/1.0.3/HLA-LA/obj \
        /opt2/hla-la/1.0.3/HLA-LA/temp \
        /opt2/hla-la/1.0.3/HLA-LA/working \
        /opt2/hla-la/1.0.3/HLA-LA/graphs \
    && cd /opt2/hla-la/1.0.3/HLA-LA/src \
    && git clone https://github.com/DiltheyLab/HLA-LA.git . \
    && git checkout fb8ae1cde7783d6cf441a80319e943f63f3b6276 \
    && ln -sf /opt2/libbamtools/lib /opt2/libbamtools/lib64 \
    && make all \
        BOOST_PATH=/opt2/libboost \
        BAMTOOLS_PATH=/opt2/libbamtools \
    && sed -i "s:^\(picard_sam2fastq_bin=\).*:\1$PICARD:" paths.ini \
    && sed -i "s:^\(workingDir=\).*:workingDir=:" paths.ini \
    && sed -i "s:^\(workingDir_HLA_ASM=\).*:workingDir_HLA_ASM=:" paths.ini \
    && /opt2/hla-la/1.0.3/HLA-LA/bin/HLA-LA --action testBinary
ENV PATH="${PATH}:/opt2/hla-la/1.0.3/HLA-LA/bin:/opt2/hla-la/1.0.3/HLA-LA/src"

# COPY Additional_B38_3.txt /opt2/hla-la/1.0.3/HLA-LA/src/additionalReferences/PRG_MHC_GRCh38_withIMGT/

################ POST #################
# Add Dockerfile and export env variables
# ADD Dockerfile /opt2/genome-seek_hla.dockerfile
# RUN chmod -R a+rX /opt2
# WORKDIR /data2