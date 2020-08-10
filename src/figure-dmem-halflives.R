#!/usr/bin/env Rscript

########################################
## filename: figure-dmem-halflives.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot figure showing half-lives in heat-treated
## dmem for all dmem options
#######################################


script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'ggplot2',    # for plotting
    'tidyr',      # for crossing()
    'cowplot',    # publication ready ggplot
    'extrafont',
    'heatinactivation'
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
material_colors <- get_params("material_colors")

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
                   "DMEM uncovered plate oven" = "uncovered\nplate oven",
                   "DMEM covered plate oven" = "covered\nplate oven",
                   "DMEM closed vial oven" = "closed vial\noven",
                   "DMEM closed vial heat block" = "closed vial\nheat block")
           )

##################################################
## calculate posterior draws for regression lines
##################################################

hl_draws <- decay_chains %>%
    spread_draws(log_half_life[experiment_id]) %>%
    add_titer_metadata(dat, "experiment_id")

hl_plot <- hl_draws %>%
    ggplot(aes(x = material,
               y = exp(log_half_life) * 60)) +
    stat_dotsinterval(quantiles = 100,
                      slab_color = "black",
                      stroke = 0,
                      fill = material_colors[["DMEM"]]) +
    scale_y_log10_mathformat() +
    theme_project() +
    theme(legend.position = "none") +
    xlab("treatment") +
    ylab("half-life (min)")


####################################
## compose full figure from panels
####################################

cat('making full figure...\n')

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          hl_plot,
          base_width = 12,
          base_height = 10)
warnings()
