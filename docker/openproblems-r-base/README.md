# openproblems-r-base Docker image

Base image: singlecellopenproblems/openproblems

OS: Debian Stretch

Python: 3.8

R: 4.0

apt packages:

* dirmngr
* ca-certificates
* gnupg
* gpgv
* gfortran
* libblas-dev
* liblapack-dev
* r-base-core=4.0

R packages:

* BiocManager
* scran
* IRKernel

Python packages:

* rpy2
* anndata2ri>=1.1
