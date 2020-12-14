if (!require(shiny, quietly=TRUE)) install.packages('shiny')
if (!require(devtools, quietly=TRUE)) install.packages('devtools')
devtools::install_version('RcppAnnoy', version = '0.0.16')
