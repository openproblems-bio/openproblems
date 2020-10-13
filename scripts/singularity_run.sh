set -x
WORKDIR=$1
SCRIPT=$2
TASKNAME=$3
FUN=$4
ARGS=${@:5}
pip install --no-cache-dir --user -q -U ${WORKDIR}/..[${TASKNAME}-${FUN}] && cd ${WORKDIR} && python3 $SCRIPT $TASKNAME $FUN ${ARGS}
