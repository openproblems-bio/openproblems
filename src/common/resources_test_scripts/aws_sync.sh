#!/bin/bash

echo "Run the command in this script manually"
exit 1

aws s3 sync "resources_test" "s3://openproblems-data" --exclude */temp_* --delete --dryrun