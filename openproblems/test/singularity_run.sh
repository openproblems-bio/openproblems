set -x
SCRIPT=$1
TASKNAME=$2
FUN=$3
ARRAY=( $@ )
LEN=${#ARRAY[@]}
ARGS=${ARRAY[@]:3:$LEN}
pip install --no-cache-dir --user -q -U ..[${TASKNAME}-${FUN}] && python3 $SCRIPT $TASKNAME $FUN ${ARGS}
