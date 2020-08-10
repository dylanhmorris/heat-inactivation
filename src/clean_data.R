#!/usr/bin/env Rscript

#####################################
## name: clean_data.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## process raw data and save cleaned
## data for use in model fitting
##
####################################


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

## functions

#' clean_data
#'
#' function to clean a raw dataset
#'
#' @param dat dataset to clean
#' 
#' @return cleaned dataset
clean_data <- function(dat,
                       min_time) {
    
        
    ## for now, do not estimate dose deposited
    ## loss rate -- maybe add if everything
    ## else works

    dat <- dat %>% filter(time >= min_time,
                          !grepl("cytotoxic", comment) &
                          virus == "SARS-CoV-2" &
                          grepl("DMEM", material))

    dat$treatment <- dat$treatment %>%
        replace_na("Untreated")

    dat <- dat %>%
        group_by(virus,
                 temperature,
                 treatment,
                 material) %>%
        mutate(experiment_id = cur_group_id())
    
    dat <- dat %>%
        arrange(material) %>%
        group_by(material) %>%
        mutate(material_id = cur_group_id())
    
    dat <- dat %>%
        group_by(experiment_id,
                 replicate,
                 time) %>%
        mutate(titer_id = cur_group_id())
    
    dat <- dat %>% arrange(titer_id)
    
    print(dat)
    return (dat)
}

#' main
#' 
#' function to execute script
#' 
#' @return nothing
main <- function(){

    ## arguments
    args <- commandArgs(trailingOnly = TRUE)
    raw_data_path <- args[1]
    out_path <- args[2]

    ## fixed params
    delim <- ";"
    min_time <- 0

    col_types <- cols_only(
        virus = col_factor(),
        material = col_factor(),
        treatment = col_character(),
        time = col_double(),
        dilution = col_integer(),
        replicate = col_integer(),
        virus_detect = col_integer(),
        temperature = col_factor(),
        humidity = col_double(),
        comment = col_character())

    cat("Reading raw data...\n")
    dat <- read_delim(raw_data_path,
                      delim = delim,
                      col_types = col_types)

    cat("cleaning raw data...\n")
    cleaned <- clean_data(dat,
                          min_time)

    cat("Writing results to ", out_path, "...\n")
    write_csv(cleaned,
              out_path)

    warnings()
}

main()
