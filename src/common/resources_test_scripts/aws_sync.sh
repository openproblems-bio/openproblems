#!/bin/bash

echo "Run the command in this script manually"
exit 1

aws s3 sync "resources_test" "s3://openproblems-data/resources_test" --exclude "*/temp*" --exclude "*/tmp*" --delete --dryrun
aws s3 sync "resources" "s3://openproblems-data/resources" --exclude */temp_* --delete --dryrun
