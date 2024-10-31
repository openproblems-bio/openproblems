#!/bin/bash

## VIASH START
par_input='_viash.yaml'
par_output='.'
## VIASH END

extra_params=( )

if [ "$par_quiet" == "true" ]; then
  extra_params+=( "--quiet" )
fi
if [ "$par_dryrun" == "true" ]; then
  extra_params+=( "--dryrun" )
fi
if [ "$par_delete" == "true" ]; then
  extra_params+=( "--delete" )
fi

if [ ! -z ${par_exclude+x} ]; then
  IFS=";"
  for var in $par_exclude; do
    unset IFS
    extra_params+=( "--exclude" "$var" )
  done
fi

function sync_s3() {
  local s3_path="$1"
  local dest_path="$2"
  AWS_EC2_METADATA_DISABLED=true \
    aws s3 sync \
    "$s3_path" \
    "$dest_path" \
    --no-sign-request \
    "${extra_params[@]}"
}

yq e \
  '.info.test_resources[] | "{type: " + (.type // "s3") + ", path: " + .path + ", dest: " + .dest + "}"' \
  "${par_input}" | \
  while read -r line; do
    type=$(echo "$line" | yq e '.type')
    path=$(echo "$line" | yq e '.path')
    dest=$(echo "$line" | yq e '.dest')

    echo "Syncing '$path' to '$dest'..."

    if [ "$type" == "s3" ]; then
      sync_s3 "$path" "$par_output/$dest"
    fi
  done
