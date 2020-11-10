#!/bin/bash
WORKDIR=$1
SCRIPT=$2
TASKNAME=$3
FUN=$4
ARRAY=( $@ )
LEN=${#ARRAY[@]}
ARGS=${ARRAY[@]:4:$LEN}
CODEDIR=$(dirname $(dirname $WORKDIR))
pip install --no-cache-dir -qq -U $CODEDIR
cd ${WORKDIR}
coverage run --parallel --source=${CODEDIR}/openproblems $SCRIPT $TASKNAME $FUN ${ARGS}
