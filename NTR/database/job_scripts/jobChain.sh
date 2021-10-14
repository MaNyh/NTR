#!/bin/bash
JOB_NUMBERS=12
DEPENDENCY=""
JOB_FILE="submit_job_hlrn.sh"

for (( c=1; c<=$JOB_NUMBERS; c++ )) ; do
    JOB_CMD="sbatch"
    if [ -n "$DEPENDENCY" ] ; then
        JOB_CMD="$JOB_CMD --dependency afterany:$DEPENDENCY"
    fi
    JOB_CMD="$JOB_CMD $JOB_FILE"
    echo -n "Running command: $JOB_CMD  "
    OUT=`$JOB_CMD`
    echo "Result: $OUT"
    DEPENDENCY=`echo $OUT | awk '{print $4}'`
done
