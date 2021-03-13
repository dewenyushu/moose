#!/bin/bash
#PBS -M dewen.yushu@inl.gov
#PBS -m be
#PBS -N sp_1
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -l place=scatter:excl
#PBS -l walltime=168:00:00
#PBS -P ldrd

# A simple script to run contact problems with different preconditioners and refinement levels

JOB_NUM=${PBS_JOBID%%\.*}

cd $PBS_O_WORKDIR

module purge
module load pbs
module load use.moose PETSc

ID=1
Tmelt=1200

DESTDIR="./output_case_${ID}"

if [ ! -d "$DESTDIR" ]
then
mkdir ${DESTDIR}
fi

INPUT="./3DP_brick.i"

OUTPUT="${DESTDIR}/case_${ID}.log"

OPTIONS="Case=${ID} T_melt=${Tmelt}"

mpiexec ~/projects/moose/modules/combined/combined-opt -i $INPUT  $OPTIONS -snes_view -log_view >${OUTPUT}
