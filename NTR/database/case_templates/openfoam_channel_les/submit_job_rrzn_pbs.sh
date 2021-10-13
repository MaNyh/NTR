#PBS -S /bin/ksh
#PBS -N <var CLUSTERJOBNAME var>
#PBS -j oe
#PBS -l nodes=<var NODES var>:ppn=<var TASKS_PER_NODE var>
#PBS -l walltime=<var WALLTIME var>
#PBS -l mem=<var MEMORY var>
#PBS -M <var JOB_EMAIL var>
#PBS -m ae
#PBS -W x=PARTITION:haku:taurus:lena:tane
. $MODULESHOME/init/ksh
module load GCC/4.9.3-2.25 OpenMPI/1.10.2 OpenFOAM/v1612+ && source $FOAM_BASH
HOST=$(paste -s -d ',' ${PBS_NODEFILE})
echo $HOST
echo $PBS_NODEFILE
cd $PBS_O_WORKDIR
mpirun --hostfile ${PBS_NODEFILE} -np <var PROCS var> pimpleFoam -parallel > log.log

