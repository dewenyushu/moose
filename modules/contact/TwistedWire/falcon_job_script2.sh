#!/bin/bash
#PBS -M dewen.yushu@inl.gov
#PBS -m be
#PBS -N Twist2
#PBS -l select=2:ncpus=36:mpiprocs=36
#PBS -l place=scatter:excl
#PBS -l walltime=10:00:00
#PBS -P neams

# A simple script to run contact problems with different preconditioners and refinement levels

JOB_NUM=${PBS_JOBID%%\.*}

cd $PBS_O_WORKDIR

module purge
module load pbs
module load use.moose PETSc/3.10.5-GCC


SOURCEDIR="./"
DESTDIR="${SOURCEDIR}output/"

if [ ! -d "$DESTDIR" ]
then
mkdir ${DESTDIR}
fi

NZ=40
Ring_Num=2
Sector=4

INPUT="${SOURCEDIR}torque_reaction_wire_genmesh.i"


OUTPUT="${DESTDIR}/nz${NZ}_ring_num${Ring_Num}_sector${Sector}.log"

mpiexec ./../contact-opt -i $INPUT nz=${NZ} ring_num=${Ring_Num} sector=${Sector} Job=2 >${OUTPUT}
