#PBS -S /bin/ksh
#PBS -N TEST_mappedDFSEM__PipeFlow
#PBS -j oe
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00,mem=16Gb
#PBS -M nyhuis@tfd.uni-hannover.de
#PBS -m ae
#PBS -W x=PARTITION:taurus:haku:lena
. $MODULESHOME/init/ksh
source ~/ModuleLoading/OpenFoam1612.sh
HOST=$(paste -s -d ',' ${PBS_NODEFILE})

cd $PBS_O_WORKDIR
rhoPimpleFoam > log.log
