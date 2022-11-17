#!/bin/bash

## VIASH START
par_input='s3://openproblems-data/resources_test'
par_output='resources_test'
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
  IFS=":"
  for var in $par_exclude; do
    unset IFS
    extra_params+=( "--exclude" "$var" )
  done
fi


# Disable the use of the Amazon EC2 instance metadata service (IMDS).
# see https://florian.ec/blog/github-actions-awscli-errors/
# or https://github.com/aws/aws-cli/issues/5234#issuecomment-705831465
export AWS_EC2_METADATA_DISABLED=true

aws s3 sync "$par_input" "$par_output" --no-sign-request "${extra_params[@]}"
