FROM ubuntu:22.04

LABEL description="MODFLOW 6 extended build with gfortran, OpenMPI, netCDF, PETSc, and a miniconda development environment."

ARG work_dir=/workdir
ARG petsc_dir=/workdir/petsc
ARG petsc_arch=arch-gcc-ompi-opt
ARG modflow_dir=/workdir/modflow6

ENV PROJECT_NAME=modflow6-extended
ENV PATH=${modflow_dir}/bin:${PATH}
ENV LD_LIBRARY_PATH=${petsc_dir}/${petsc_arch}/lib:/usr/local/lib:/usr/lib:${LD_LIBRARY_PATH}
ENV MPI_VERSION=5.0.7
ENV MPI_SHORT_VER=v5.0
ENV PETSC_VERSION=3.22.4

# install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    bash \
    gcc \
    gfortran \
    g++ \
    make \
    python3 \
    zlib1g-dev \
    pkg-config \
    build-essential \
    libnetcdf-dev \
    libnetcdff-dev \
    netcdf-bin
	
# install miniconda
WORKDIR ${work_dir}
RUN mkdir -p /opt/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda3/miniconda.sh
RUN bash /opt/miniconda3/miniconda.sh -b -u -p /opt/miniconda3
RUN rm /opt/miniconda3/miniconda.sh
ENV PATH="/opt/miniconda3/bin:$PATH"

# configure netCDF
RUN nc-config --all
RUN nf-config --all

# download, configure and build OpenMPI
WORKDIR ${work_dir}
RUN wget https://download.open-mpi.org/release/open-mpi/${MPI_SHORT_VER}/openmpi-${MPI_VERSION}.tar.gz
RUN gunzip -c openmpi-${MPI_VERSION}.tar.gz | tar xf -
WORKDIR ${work_dir}/openmpi-${MPI_VERSION}
RUN ./configure --prefix=/usr/local
RUN make all install
ENV LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}
ENV PATH=/usr/local/bin:${PATH}

# download, configure and build PETSc
WORKDIR ${work_dir}
RUN wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz
RUN gunzip -c petsc-${PETSC_VERSION}.tar.gz | tar xf -
RUN mv petsc-${PETSC_VERSION} ${petsc_dir}
ENV PETSC_DIR=${petsc_dir}
ENV PETSC_ARCH=${petsc_arch}
WORKDIR ${petsc_dir}
RUN ./configure --download-fblaslapack --with-debugging=0 --with-mpi-dir=/usr/local COPTFLAGS=-O2 CXXOPTFLAGS=-O2 FOPTFLAGS=-O2
RUN make all
ENV LD_LIBRARY_PATH=${PETSC_DIR}/${PETSC_ARCH}/lib:${LD_LIBRARY_PATH}
ENV PKG_CONFIG_PATH=${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:${PKG_CONFIG_PATH}

# checkout, build and install MODFLOW 6
WORKDIR ${work_dir}
RUN git clone https://github.com/MODFLOW-ORG/modflow6
WORKDIR ${work_dir}/modflow6
# TODO: how to avoid the conda boilerplate in these commands?
# no suggestions in https://stackoverflow.com/a/62674910...
RUN conda init bash && . ~/.bashrc && conda env create -f environment.yml
RUN conda init bash && \
    . ~/.bashrc && \
    conda activate modflow6 && \
    meson setup builddir --prefix=$(pwd) --libdir=bin -Ddebug=false -Dextended=true && \
    meson install -C builddir
