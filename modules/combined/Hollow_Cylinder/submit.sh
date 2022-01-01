#!/bin/bash
#PBS -M dewen.yushu@inl.gov
#PBS -m be
#PBS -N cylinder
#PBS -l select=1:ncpus=48:mpiprocs=32
#PBS -l place=scatter:excl
#PBS -l walltime=168:00:00
#PBS -P ne_ldrd

# A simple script to run contact problems with different preconditioners and refinement levels

JOB_NUM=${PBS_JOBID%%\.*}

cd $PBS_O_WORKDIR

module purge
module load pbs
module load use.moose PETSc

# DESTDIR="./log_files"
#
# if [ ! -d "$DESTDIR" ]
# then
# mkdir ${DESTDIR}
# fi

INPUT="master_app_mechanical.i"

# OUTPUT="${DESTDIR}/log.txt"

OPTIONS=" "

# mpiexec ~/projects/moose_AM/modules/combined/combined-opt -i $INPUT $OPTIONS -log_view >${OUTPUT}
mpiexec ~/projects/moose_AM/modules/combined/combined-opt -i $INPUT $OPTIONS
