#!/bin/bash

##########
#
#   PBS job submission information
#
##########

# Choose nodes (or number of nodes) on which to run
#####PBS -l nodes=node22:ppn=3+node23:ppn=3+node24:ppn=3
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00

# priority
#PBS -q immediate

# Standard out, error
#PBS -o mx.out
#PBS -e mx.err

# Export all environment variables to the job
#PBS -V


##########
#
#   If on PBS, get to correct directory
#
##########
if test -n "$PBS_NODEFILE"; then
  cat $PBS_NODEFILE
  export STARTDIR=$PWD
  cd ${PBS_O_WORKDIR}
  echo STARTDIR = $STARTDIR
  echo PBS_O_WORKDIR = $PBS_O_WORKDIR
fi



##########
#
#   Run!
#
##########

MXBIN=../build/src/maxwell
MXIN=phc-sapph-r0.37-036.mx

mpiexec ${MXBIN} --infile=${MXIN}
