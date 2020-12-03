#!/bin/bash
set -x
WORKDIR=$1
SCRIPT=$2
TASKNAME=$3
FUN=$4
ARRAY=( $@ )
LEN=${#ARRAY[@]}
ARGS=${ARRAY[@]:4:$LEN}
CODEDIR=$(dirname $(dirname $WORKDIR))
python3 -m pip install --no-cache-dir -qq -U $CODEDIR
python3 -m pip install -qq -U coverage
cd ${WORKDIR}
python3 -m coverage run --parallel --source=${CODEDIR}/openproblems $SCRIPT $TASKNAME $FUN ${ARGS}
