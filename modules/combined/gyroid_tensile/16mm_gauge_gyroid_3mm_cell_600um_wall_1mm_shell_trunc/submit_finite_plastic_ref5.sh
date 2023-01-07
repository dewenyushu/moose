#!/bin/bash
#PBS -M dewen.yushu@inl.gov
#PBS -m be
#PBS -N trunc_ref5
#PBS -l select=2:ncpus=48:mpiprocs=8
#PBS -l place=scatter:excl
#PBS -l walltime=40:00:00
#PBS -P ne_ldrd

JOB_NUM=${PBS_JOBID%%\.*}

cd $PBS_O_WORKDIR

module purge
module load use.moose PETSc

SOURCEDIR="./"
DESTDIR="./output/"
#
# if [ ! -d "$DESTDIR" ]
# then
# mkdir ${DESTDIR}
# fi

INPUT="${SOURCEDIR}input_disp_control_plastic_ref5.i"

OUTPUT="${DESTDIR}input_disp_control_plastic_ref5.log"

mpiexec ~/projects/light_weight_ldrd_moose/modules/combined/combined-opt -i $INPUT --recover -snes_view -log_view >${OUTPUT}
