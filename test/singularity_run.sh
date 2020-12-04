#!/bin/bash
WORKDIR=$1
SCRIPT=$2
TASKNAME=$3
FUN=$4
ARRAY=( $@ )
LEN=${#ARRAY[@]}
ARGS=${ARRAY[@]:4:$LEN}
CODEDIR=$(dirname $WORKDIR)
export PYTHONPATH=$CODEDIR

python3 -m pip install --upgrade "pip<=21.0"
python3 -m pip install --use-deprecated=legacy-resolver --upgrade-strategy=only-if-needed --no-cache-dir -U $CODEDIR
python3 -m pip install --use-deprecated=legacy-resolver -U coverage

cd ${WORKDIR}
python3 -m coverage run --source=${CODEDIR}/openproblems $SCRIPT $TASKNAME $FUN ${ARGS}
