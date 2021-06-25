#PBS -S /bin/ksh 
#PBS -N TEST_mappedDFSEM__PipeFlow 
#PBS -j oe 
#PBS -l nodes=1:ppn=12
#PBS -l walltime=10:00:00,mem=16Gb 
#PBS -M nyhuis@tfd.uni-hannover.de 
#PBS -m ae 
#PBS -W x=PARTITION:taurus:haku:lena
. $MODULESHOME/init/ksh 
source ~/ModuleLoading/OpenFoam1612.sh 
HOST=$(paste -s -d ',' ${PBS_NODEFILE}) 
echo $HOST 
echo $PBS_NODEFILE

cd $PBS_O_WORKDIR 
mpirun --hostfile ${PBS_NODEFILE} -np 12 rhoPimpleFoam -parallel > ${PBS_JOBNAME}.${PBS_JOBID}.out
#rhoPimpleFoam > ssp.${PBS_JOBNAME}.${PBS_JOBID}.out
