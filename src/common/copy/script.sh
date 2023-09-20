#!/bin/bash

set -eo pipefail

## VIASH START
par_input="input.txt"
par_output="output.txt"
## VIASH END

parent=`dirname "$par_output"`
if [[ ! -d "$parent" ]]; then
  mkdir -p "$parent"
fi

cp -r "$par_input" "$par_output"