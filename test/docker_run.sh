#!/bin/bash
WORKDIR=$1
SCRIPT=$2
ARRAY=( "$@" )
LEN=${#ARRAY[@]}
ARGS=( "${ARRAY[@]:2:$LEN}" )
CODEDIR=$(dirname "$WORKDIR")
export PYTHONPATH="$WORKDIR"

if [ ! -f ~/.install_complete ]; then
  python3 -m pip install --upgrade "pip<=21.0"
  python3 -m pip install --use-deprecated=legacy-resolver --upgrade-strategy=only-if-needed --no-cache-dir --editable "$CODEDIR"
  python3 -m pip install --use-deprecated=legacy-resolver -U coverage
  FREEZE="$(python3 -m pip freeze)"
  if echo "$FREEZE" | grep -q annoy; then
    python3 -m pip install --force "$(echo "$FREEZE" | grep annoy)"
  fi
  touch ~/.install_complete
fi

cd "${CODEDIR}" || exit
python3 -m coverage run --parallel --source=openproblems "${SCRIPT}" "${ARGS[@]}"
