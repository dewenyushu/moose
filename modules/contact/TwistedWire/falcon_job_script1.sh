#!/bin/bash
#PBS -M dewen.yushu@inl.gov
#PBS -m be
#PBS -N Twist1
#PBS -l select=1:ncpus=36:mpiprocs=36
#PBS -l place=scatter:excl
#PBS -l walltime=50:00:00
#PBS -P neams

# A simple script to run contact problems with different preconditioners and refinement levels

JOB_NUM=${PBS_JOBID%%\.*}

cd $PBS_O_WORKDIR

module purge
module load pbs
module load use.moose PETSc/3.10.5-GCC


SOURCEDIR="./"
DESTDIR="${SOURCEDIR}stroing_scaling_output/36cores/"

if [ ! -d "$DESTDIR" ]
then
mkdir ${DESTDIR}
fi

INPUT="${SOURCEDIR}torque_reaction_wire_genmesh.i"

OUTPUT="${DESTDIR}/36cores.log"

mpiexec ./../contact-opt -i $INPUT Job=1 >${OUTPUT}
