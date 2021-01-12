# openproblems-r-base Docker image

Base image: singlecellopenproblems/openproblems

OS: Debian Stretch

Python: 3.7

R: 3.6

apt packages:

* dirmngr
* ca-certificates
* gnupg
* gpgv
* gfortran
* libblas-dev
* liblapack-dev
* r-base-core=3.6

R packages:

* BiocManager
* scran

Python packages:

* rpy2<3.4.0
* scIB
* anndata2ri
