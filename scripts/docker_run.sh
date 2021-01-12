set -x
WORKDIR=$1
SCRIPT=$2
ARGS=${@:3}
cd ${WORKDIR}
python3 $SCRIPT ${ARGS}
