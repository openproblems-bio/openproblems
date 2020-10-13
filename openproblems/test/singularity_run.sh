#!/bin/bash
WORKDIR=$1
SCRIPT=$2
TASKNAME=$3
FUN=$4
ARRAY=( $@ )
LEN=${#ARRAY[@]}
ARGS=${ARRAY[@]:4:$LEN}
pip install --no-cache-dir -qq -U $(dirname $(dirname $WORKDIR))[${TASKNAME}-${FUN}] && cd ${WORKDIR} && python3 $SCRIPT $TASKNAME $FUN ${ARGS}
