## VIASH START
par_input=s3://openproblems-nextflow/work/f6/8565066aee4771cc2790b92b4ac660
par_aws_profile=op
par_aws_credentials=~/.aws/credentials
par_output=debug
## VIASH END

if [ -n "$par_aws_credentials" ]; then
  export AWS_SHARED_CREDENTIALS_FILE="$par_aws_credentials"
fi

if [ -n "$par_aws_profile" ]; then
  export AWS_PROFILE="$par_aws_profile"
fi

aws s3 sync "$par_input" "$par_output"

cd "$par_output"

# Replace the miniconda aws with the system aws
sed -i 's#/home/ec2-user/miniconda/bin/aws#aws#g' .command.run
# sed -i 's# s3 cp # s3 sync #g' .command.run

bash .command.run nxf_stage
