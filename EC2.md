# AWS EC2 Usage Instructions

The following instructions give a step-by-step guide to launching an AWS virtual machine
with all the required prerequisites to run `openproblems`.

## Code of conduct

**Please be respectful of our finite resources!**

* The use of the `openproblems` AWS account is a privilege, not a right.
* This privilege is given solely for the purposes of testing methods with
  `openproblems-cli test`.
* Developers who have their own compute resources should use them; please help us
  conserve our resources for those who need them.
* If developers are found to be using resources irresponsibly, we may have to revoke
  this privilege.

## Requirements

* a Unix-based OS (Mac or Linux), though you should be
able to amend the commands for use on Windows (or consider [Windows Subsystem for
Linux](https://docs.microsoft.com/en-us/windows/wsl/install)).
* The [AWS CLI](https://aws.amazon.com/cli/)
* [jq](https://stedolan.github.io/jq/download/)

## Instructions

The following instructions are for `bash`, other shell users may need to modify commands
slightly.

First, if you have recieved openproblems AWS credentials, configure AWS to use them.
Note: `openproblems` uses `us-west-2` as default region. If you have other AWS accounts,
you can configure AWS with multiple accounts by using the `AWS_PROFILE` environment
variable.

```shell
export AWS_PROFILE=openproblems
aws configure
```

Second, create a key pair (only do this once):

```shell
KEY_NAME="my_openproblems_key" # name this whatever you like, but it must be unique
aws ec2 create-key-pair --key-name $KEY_NAME --key-format pem \
--query "KeyMaterial" --output text > ${KEY_NAME}.pem
chmod 400 ${KEY_NAME}.pem
```

Now, create an instance with your key pair:

```shell
OWNER_NAME="this_is_your_name"
AWS_EC2_INSTANCE_TYPE="t2.micro"
INSTANCE_ID=$(
aws ec2 run-instances --count 1 --image-id ami-01219569b1bbf9fb2 \
  --instance-type $AWS_EC2_INSTANCE_TYPE --key-name $KEY_NAME \
  --security-group-ids sg-002d2b9db29bb43dd \
  --tag-specifications "ResourceType=instance,Tags=[{Key=owner,Value=${OWNER_NAME}}]" |
  jq '.["Instances"][0]["InstanceId"]' |
  tr -d '"'
)
```

Get the public DNS address for your instance

```shell
sleep 30 # wait for boot
PUBLIC_DNS_NAME=$(
aws ec2 describe-instances --instance-id $INSTANCE_ID |
  jq '.["Reservations"][0]["Instances"][0]["PublicDnsName"]' |
  tr -d '"'
)
```

Now you can SSH into your instance:

```shell
# check the status of your instance
aws ec2 describe-instance-status --instance-id ${INSTANCE_ID} | \
jq '.["InstanceStatuses"][0]["SystemStatus"]'
ssh -i ${KEY_NAME}.pem ubuntu@${PUBLIC_DNS_NAME}
```

The instance will by default contain all dependencies to use `openproblems`. You can
run `openproblems` with

```shell
git clone https://github.com/openproblems-bio/openproblems
cd openproblems
sudo docker run \
  -v $(pwd):/usr/src/singlecellopenproblems -v /tmp:/tmp \
  -it singlecellopenproblems/openproblems bash
openproblems-cli --help
```

For more information on using the CLI, see
[CONTRIBUTING.md](CONTRIBUTING.md#testing-method-performance).

When you are done, make sure to shut down your instance:

```shell
aws ec2 terminate-instances --instance-ids ${INSTANCE_ID}
```

Finally, make sure you don't have any instances left running:

```shell
aws ec2 describe-instances --filters "Name=tag:owner,Values=${OWNER_NAME}"
```
