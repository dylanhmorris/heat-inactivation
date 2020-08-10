#!/usr/bin/env Rscript

########################################
## filename: table-dmem-halflives.R
## author: Dylan Morris <dhmorris@princeton.edu>
## create table of half-life posterior medians
## and 95% CI
#######################################

script_packages <- c(
    'heatinactivation', # project functions
    'rstan',      # stan interface
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'xtable'      # to create LaTeX table
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}

#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly=TRUE)
decay_data_path <- args[1]
decay_results_path <- args[2]
titer_path <- args[3]
outpath <- args[4]

decay_chains <- readRDS(decay_results_path)
dat <- read_csv(decay_data_path,
                col_types = cols())

#################################
## overall plot styling
#################################

material_order <- c(
    "DMEM uncovered plate oven",
    "DMEM covered plate oven",
    "DMEM closed vial oven",
    "DMEM closed vial heat block",
    "DMEM")

dat$material <- factor(
    dat$material,
    levels = material_order)

dat <- dat %>%
    mutate(material = material %>%
               recode_factor(
                   "DMEM uncovered plate oven" = "Uncovered\nplate oven",
                   "DMEM covered plate oven" = "Covered\nplate oven",
                   "DMEM closed vial oven" = "Closed vial\noven",
                   "DMEM closed vial heat block" = "Closed vial\nheat block")
           )

##################################################
## calculate half-life table
##################################################

hl_draws <- decay_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(dat, "experiment_id")

hl_table <- hl_draws %>%
    mutate(hl_min = exp(log_half_life) * 60) %>%
    group_by(material) %>%
    summarise(median = median(hl_min),
              q025 = quantile(hl_min, 0.025),
              q975 = quantile(hl_min, 0.975)) %>%
    xtable()


####################################
## save table
####################################

cat('saving table to', outpath, '...\n')

print(hl_table, file = outpath)

warnings()
