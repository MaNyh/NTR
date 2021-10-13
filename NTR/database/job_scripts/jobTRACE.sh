#!/bin/bash
###############################################################################
###                              job_TRACE.sh                               ###
###			start vom run-directory auf Knoten f?r openmpi	    ###
###                      28 April 2009 Florian Herbst                       ###
###                                                                         ###
###############################################################################

JOBNAME=<var JOB_NAME var>
EMAIL=<var JOB_MAIL var>
NODES=<var JOB_NODES var>
#NODES="knoten-10.tfd.uni-hannover.de"
NPROCS=<var JOB_PPN var>


################################################################################
#DON'T CHANGE PART BELOW
################################################################################
WALLTIME=<var JOB_RUNTIME var>
QUEUE="all"
BASENAME="TRACE.cgns"
CONTROLFILE="TRACE_control.input"
#
let CPUS=($NPROCS \* $NODES)

#CPUS=$NPROCS
SOFTWARE_DIR="/sw/TRACE"
POST_DIR="$SOFTWARE_DIR/post45"
TRACE_EXECUTABLE="$SOFTWARE_DIR/trace_7.1.22_openmpi_MTU/TRACE"
BALANCE_FILE="BALANCE_"$CPUS"PROC"
MPIDIR="$SOFTWARE_DIR/openmpi-1.3.1"
RUNDIR=$PWD
echo $CPUS
#
cd ../
WORKDIR=$PWD
cd input/
$POST_DIR/post45 -cgns $BASENAME -pp $CPUS
cd $RUNDIR
#
#set TRACE environment variables
#TRACE3D_PARSE_FILE_IN="$WORKDIR/input/TRACE_control.input"
#TRACE3D_INFO_FILE_IN="$WORKDIR/input/$BALANCE_FILE"
#################################################################################
#Write tempscript to start TRACE
rm -f $RUNDIR/tempscript.sh
touch $RUNDIR/tempscript.sh
chmod +x $RUNDIR/tempscript.sh
echo -e "#!/bin/bash" >> $RUNDIR/tempscript.sh;
echo -e "cd $RUNDIR" >> $RUNDIR/tempscript.sh;
#Mit Controlfileecho -e "$MPIDIR/bin/mpiexec $TRACE_EXECUTABLE -cntl $WORKDIR/input/$CONTROLFILE -cgns $WORKDIR/input/$BASENAME -np $CPUS -lb $WORKDIR/input/$BALANCE_FILE" >> $RUNDIR/tempscript.sh;
echo -e "$MPIDIR/bin/mpiexec $TRACE_EXECUTABLE -cgns $WORKDIR/input/$BASENAME -lb $WORKDIR/input/$BALANCE_FILE -o TRACE.lst. -pc" >> $RUNDIR/tempscript.sh;
#################################################################################
#submit job to queue
qsub -N $JOBNAME -l nodes=$NODES:ppn=$NPROCS,walltime=$WALLTIME -q $QUEUE -M $EMAIL -m abe -V $RUNDIR/tempscript.sh

