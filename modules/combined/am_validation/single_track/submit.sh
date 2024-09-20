#!/bin/bash
#PBS -m abe
#PBS -M dewen.yushu@inl.gov
#PBS -N val_single
#PBS -l select=1:ncpus=48:mpiprocs=48
#PBS -l place=scatter:excl
#PBS -l walltime=100:00:00
#PBS -P moose

JOB_NUM=${PBS_JOBID%%\.*}

cd $PBS_O_WORKDIR

module purge
source activate base
conda activate moose

mpiexec /home/yushdewe/projects/moose-am/modules/combined/combined-opt -i sample.i

