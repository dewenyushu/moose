#!/bin/bash
#PBS -m abe
#PBS -M dewen.yushu@inl.gov
#PBS -N case_3
#PBS -l select=3:ncpus=48:mpiprocs=48
#PBS -l place=scatter:excl
#PBS -l walltime=100:00:00
#PBS -P moose

JOB_NUM=${PBS_JOBID%%\.*}

cd $PBS_O_WORKDIR

module purge
source activate base
conda activate moose

mpiexec /home/yushdewe/projects/moose-am/modules/combined/combined-opt -i thermal.i

