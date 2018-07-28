# setup a virtual docker environment for building PFLOTRAN

FROM ubuntu:latest
RUN apt-get update -q
RUN apt-get install -y python
RUN apt-get install -y gfortran
RUN apt-get install -y g++
RUN apt-get install -y git
RUN apt-get install -y make
RUN apt-get install -y cmake
RUN apt-get install -y valgrind

RUN mkdir /app
WORKDIR /app

ADD . /app

ENV PETSC_DIR /app/petsc
ENV PETSC_ARCH gnu-c-debug
RUN git clone https://bitbucket.org/petsc/petsc && cd petsc && git checkout xsdk-0.2.0 && ./config/configure.py --PETSC_ARCH=gnu-c-debug --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --CFLAGS='-g -O0' --CXXFLAGS='-g -O0' --FFLAGS='-g -O0 -Wno-unused-function' --with-clanguage=c --with-shared-libraries=0 --with-debugging=1 --download-mpich=yes --download-hdf5=yes --with-valgrind=1 --download-parmetis=yes --download-metis=yes --download-fblaslapack=yes --with-c2html=0 && make all

RUN cd /app/src/pflotran && make pflotran && make test


