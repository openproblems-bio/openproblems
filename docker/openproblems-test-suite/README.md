# openproblems-test-suite Docker image

Base image: singlecellopenproblems/openproblems-r-base

OS: Debian Stretch

Python: 3.7

R: 3.6

Apt packages:
* docker-ce
* singularity-container

This docker image is used exclusively for running the openproblems test suite. Run as follows:

```
sudo docker run --privileged=true -d --mount type=bind,source=$(pwd),target=/opt/openproblems -it singlecellopenproblems/openproblems-test-suite /sbin/init
# cd /opt/openproblems && service docker start && python3.7 -m nose2 openproblems.test
```
