#!/usr/bin/env bash
# This script is used to set up our custom Ubuntu AMI for running `test_benchmark` on an AWS self-hosted runner

# Setup
set -exo pipefail
sudo apt-get update
sudo apt-get install -y libhdf5-dev pandoc gfortran libblas-dev liblapack-dev libedit-dev llvm-dev unzip curl

# Install Python
sudo apt-get install -y python3.8
curl https://bootstrap.pypa.io/get-pip.py > get-pip.py
sudo python3 get-pip.py
sudo pip3 install --upgrade pip

# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# Install Nextflow
sudo apt-get install -y openjdk-17-jdk
mkdir -p /tmp/nextflow
cd /tmp/nextflow || exit
wget -qO- get.nextflow.io | bash
sudo cp /tmp/nextflow/nextflow /usr/local/bin/nextflow
sudo chmod 755 /usr/local/bin/nextflow
cd || exit
nextflow -version

# Install AWS & S3FS
sudo apt-get install -y s3fs
mkdir /tmp/awscli
cd /tmp/awscli || exit
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip -q awscliv2.zip
sudo ./aws/install || sudo ./aws/install --update
aws --version

# Install Python dependencies
sudo pip3 install -U wheel setuptools
sudo pip3 install --ignore-installed --upgrade PyYAML
sudo pip3 install "openproblems[evaluate] @ git+https://github.com/openproblems-bio/openproblems"
sudo pip3 uninstall -y openproblems
