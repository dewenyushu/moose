#!/bin/bash
#PBS -M dewen.yushu@inl.gov
#PBS -m be
#PBS -N grip_load
#PBS -l select=2:ncpus=48:mpiprocs=96
#PBS -l place=scatter:excl
#PBS -l walltime=100:00:00
#PBS -P ne_ldrd

module purge
module load mpich

JOB_NUM=${PBS_JOBID%%\.*}

export GMSH_DIR=${HOME}/projects/gmsh/install/usr/local
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GMSH_DIR}/lib64
export SINGULARITYENV_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/lib64
export SINGULARITY_BIND="/etc/profile.d/modules.sh:/etc/profile.d/modules.sh:ro,/apps:/apps:ro,/hpc-common:/hpc-common:ro"

cd $PBS_O_WORKDIR
MPI_OPTIONS=""
mpiexec ${MPI_OPTIONS} singularity exec /apps/containers/rocky_ib.sif ${HOME}/projects/moose/modules/combined/combined-opt -i train.i
