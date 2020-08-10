#!/usr/bin/env Rscript

########################################
## filename: figure-predictive-check.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot predictive checks for heat inactivation
#######################################

suppressPackageStartupMessages(library(heatinactivation))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(extrafont))


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)

if( length(args) > 3) {
    data_path <- args[1]
    check_results_path <- args[2]
    titers_path <- args[3]
    outpath <- args[4]
} else if (length(args) == 3) {
    data_path <- args[1]
    check_results_path <- args[2]
    titers_path <- args[2]
    outpath <- args[3]
} else {

    stop("Must at least supply data path, check results path, and output path")
}

## read data / style files

cat("reading data (this may take a while)...\n")
check_chains <- readRDS(check_results_path)
titer_ests_chains <- readRDS(titers_path)

dat <- read_csv(data_path,
                col_types = cols())

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
                   "DMEM uncovered plate oven" = "Uncovered plate oven",
               "DMEM covered plate oven" = "Covered plate oven",
               "DMEM closed vial oven" = "Closed vial oven",
               "DMEM closed vial heat block" = "Closed vial heat block")
           )


cat("data read succesfully!\n")


#################################
## overall plot styling
#################################
set.seed(34634) # reproducible! (since we use random draws)

gen_params <- get_params("general_param_list")
general_palette <- get_params("general_palette")
material_colors <- get_params("material_colors")

text_size = 40
detection_linesize = 0.75
titer_ylab <- expression("virus titer (TCID"[50] * "/mL media)")
xlim <- c(0, max(dat$time))
ylim <- c(1e-1, 1e6)
LOD_log10_per_ml <- 0.5
LOD <- 10^LOD_log10_per_ml
line_alpha <- gen_params[["line_alpha"]]
interval_grey <- general_palette[["interval_grey"]]

conversion_factor <- 1
n_lines <- 10


##################################################
## calculate posterior draws for regression lines
##################################################

cat("extracting draws for decay rates / intercepts (this may also take a while)...\n")

tidy_draws <- check_chains %>%
    spread_draws(sampled_titers_pred[titer_id])


cat("extracting positive wells...\n")
pos_wells <- dat %>%
    group_by(titer_id) %>%
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))

print(pos_wells)

cat("extracting decay rates...\n")
tidy_draws <- tidy_draws %>%
    add_titer_metadata(dat, "titer_id")

cat('extracting titer estimates...\n')
titer_ests_draws <- titer_ests_chains %>%
    spread_draws(log10_titer[titer_id]) %>%
    add_titer_metadata(dat, "titer_id") %>%
    inner_join(pos_wells,
               by = "titer_id") %>%
    mutate(detectable = n_pos > 1)


cat("converting units and setting undetectables to detection limit\n")

## convert from TCID50/(0.1mL) to TCID50/mL
## and visualize 0 positive well titers at
## the traditional LOD
titer_ests_draws <- titer_ests_draws %>%
    mutate(log10_titer_per_ml = ifelse(
               detectable,
               log10_titer + 1,
               LOD_log10_per_ml)) %>%
    arrange(desc(time), desc(material)) ## so masks are on top 


## adding one to convert to per mL from per 0.1 mL

shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))

scale <- 2 * sqrt(xlim[2])
if(grepl("posterior", outpath))
    scale <- xlim[2]

panel <- tidy_draws %>%
    ggplot(aes(x = time,
               y = 10^(1 + sampled_titers_pred),
               fill = material)) +
    geom_hline(aes(yintercept = LOD),
               size = 2,
               linetype = "dotted") +
    stat_eye(color = "black",
             interval_alpha = 0,
             point_alpha = 0,
             scale = scale,
             size = 15,
             alpha = 0.5) +
    stat_pointinterval(
        mapping = aes(x = time,
                      y = 10^log10_titer_per_ml,
                      shape = detectable,
                      fill = material,
                      group = titer_id),
        size = 4,
        fatten_point = 3,
        data = titer_ests_draws) +
    scale_fill_manual(values = unlist(material_colors)) +
    scale_fill_manual(values = unlist(material_colors),
                       aesthetics = "point_fill") +
    scale_color_manual(values = unlist(material_colors)) +
    shape_scale + 
    scale_y_log10_mathformat(expand = c(0, 0)) +
    coord_cartesian(ylim = ylim,
                    xlim = xlim) +
    facet_wrap(vars(material))

## styling: no facet labels because is background plot
panel <- panel +
    theme_project(base_size = text_size) +
    xlab("time (hrs)") +
    ylab(titer_ylab) +
    theme(legend.position = "none")

####################################
## compose full figure from panels
####################################

labeled_panel_theme <- theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = text_size),
    strip.placement = "outside",
    strip.switch.pad.grid = unit("0.5", "in"),
    plot.subtitle = element_text(hjust = 0.5),
    plot.tag=element_text(angle=-90,
                          size = text_size),
    plot.tag.position=c(1.05, 0.5))

left_margin <- theme(
    plot.margin = margin(b = 3, t = 1, l = 1, r = 0, unit = "cm"))
    
cat('making full figure...\n')

full_fig <- panel + labeled_panel_theme + left_margin

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          full_fig,
          base_height = 15,
          base_asp = 1.2)
warnings()
