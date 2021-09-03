#!/bin/bash -l
#SBATCH --job-name=fLE_highBD_rotating2012_roundTE
#SBATCH --nodes=<var NODES var>
#SBATCH --ntasks-per-node=<var TASKS_PER_NODE var>
#SBATCH --mem-per-cpu=2G
#SBATCH --time=100:00:00

#SBATCH --mail-user=<var HLRN_JOB_EMAIL var>
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output batch.out
#SBATCH --error batch.err

. $MODULESHOME/init/ksh
source ~/ModuleLoading/OpenFoam1612.sh
HOST=$(paste -s -d ',' ${PBS_NODEFILE})
echo $HOST
echo $PBS_NODEFILE
cd $PBS_O_WORKDIR

mpirun -np <var PROCS var> rhoPimpleFoam -case ./ -parallel >& log.log
