#!/usr/bin/env Rscript

#####################################
## name: install_needed_packages.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## installs needed packages for
## reproducing SARS-CoV-2 heat
## inactivation study
##
####################################

install_if_absent <- function(package_name){
    if (!(package_name %in% installed.packages()))
        install.packages(pkgs = package_name,
                         repos = "http://cloud.r-project.org")
    else
        cat(sprintf("Package %s already installed\n", package_name))
}

install_local_if_absent <- function(package_name){
    if (!(package_name %in% installed.packages())){
        
        cat(sprintf("Installing local package %s...\n",
                    package_name))
        devtools::install(package_name)

    } else {

        cat(sprintf("Package %s already installed\n", package_name))

    }
}


## install devtools from CRAN
install_if_absent("devtools")

## install custom package for this project
## and any missing dependencies
library(devtools)
install_local_if_absent("heatinactivation")


## install fonts

r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos = r)

library(extrafont)

font_install('fontcm',
             prompt = FALSE)

