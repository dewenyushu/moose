#!/bin/bash
#PBS -m abe
#PBS -M dewen.yushu@inl.gov
#PBS -N case_4
#PBS -l select=10:ncpus=48:mpiprocs=48
#PBS -l walltime=100:00:00
#PBS -P moose

cd $PBS_O_WORKDIR
module purge
conda init
conda activate moose

mpirun ~/projects/moose-am/modules/combined-opt -i thermal.i

