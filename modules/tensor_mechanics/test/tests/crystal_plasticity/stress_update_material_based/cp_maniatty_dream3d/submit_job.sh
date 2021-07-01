#!/bin/bash
#PBS -M dewen.yushu@inl.gov
#PBS -m abe
#PBS -N cp_149xstal
#PBS -l select=2:ncpus=48:mpiprocs=48
#PBS -l walltime=24:00:00
#PBS -P neams

cd $PBS_O_WORKDIR

module purge
module load pbs
module load use.moose PETSc

mpiexec $HOME/projects/moose/modules/tensor_mechanics/tensor_mechanics-opt -i cp_mat_based.i -log_view >> cp_mat_based.log
