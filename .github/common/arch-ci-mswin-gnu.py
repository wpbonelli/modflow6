#!/usr/bin/env python3

configure_options = [
  '--with-mpiexec=/C/Program Files/Microsoft MPI/Bin/mpiexec',
  '--with-shared-libraries=0',
  # '--with-batch=1',
  '--with-debugging=0',
  # not using -g so that the binaries are smaller
  'COPTFLAGS=-O',
  'FOPTFLAGS=-O',
  'CXXOPTFLAGS=-O',
  '--with-visibility=0',
  'FFLAGS=-fno-backtrace -ffree-line-length-0',
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
