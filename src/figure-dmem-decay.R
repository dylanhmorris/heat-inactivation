#!/usr/bin/env Rscript

########################################
## filename: figure-dmem-decay.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot figure showing decay in heat-treated
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
titer_chains <- readRDS(titer_path)
dat <- read_csv(decay_data_path,
                col_types = cols())

#################################
## overall plot styling
#################################
set.seed(989327) # reproducible! (since we use random draws)
n_lines <- 10
line_alpha <- 0.1
material_colors <- get_params("material_colors")
titer_ylab <- expression("virus titer (TCID"[50] * "/mL media)")

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

##################################################
## calculate posterior draws for regression lines
##################################################

plot_times <- tibble(time = seq(0, 2, length.out = 1000))

## get needed draws and add human readable names
int_draws <- decay_chains %>%
    spread_draws(intercept[titer_id]) %>%
    add_titer_metadata(dat)

decay_draws <- decay_chains %>%
    spread_draws(decay_rate[experiment_id])

draws <- int_draws %>%
    inner_join(decay_draws,
               by = c(".draw", "experiment_id"))

pos_wells <- dat %>%
    group_by(titer_id) %>%
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))    


titer_draws <- titer_chains %>%
    spread_draws(log10_titer[titer_id]) %>%
    add_titer_metadata(dat) %>%
    inner_join(pos_wells,
               by = "titer_id") %>%
    mutate(detectable = n_pos > 1)


## convert from TCID50/(0.1mL) to TCID50/mL
## and visualize 0 positive well titers at
## the traditional LOD
LOD_log10_per_ml = 0.5
LOD = 10^LOD_log10_per_ml
titer_draws <- titer_draws %>%
    mutate(log10_titer_per_ml = ifelse(
               detectable,
               log10_titer + 1,
               LOD_log10_per_ml)) %>%
    ungroup() %>%
    arrange(desc(time), desc(material))


###################################
## plot panel showing raw surface
## data
###################################
cat('plotting raw data...\n')

###################################
## plot panel showing fit of
## regression lines to real data
###################################
cat('plotting regression lines...\n')
## draw n_lines random regression lines
func_samples <- draws %>%
    group_by(titer_id) %>%
    sample_n(n_lines) %>%
    ungroup()

## annotate lines so that each
## has a unique id for ggplot overplotting
## (else two lines from the same draw but
## different replicates can get confused
## with each other)
func_samples <- func_samples %>%
    mutate(line_id = as.numeric(rownames(func_samples)))

## cross product decay_rates with x (time) values
## and calculate y (titer) values
cat('setting up x values...\n')
to_plot <- func_samples %>%
    rename(t = time) %>%
    crossing(plot_times)

## adding one to convert to per mL from per mL/10
to_plot <- to_plot %>%
    mutate(predicted_titer = 10^(1 + intercept - decay_rate * time))


shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))

fit_panel <- to_plot %>%
    ggplot(aes(x = time,
               y = predicted_titer,
               group = line_id)) +
    geom_hline(aes(yintercept = LOD),
               size = 2,
               linetype = "dotted") +
    geom_line(
        aes(color = material),
        alpha = line_alpha) +
    stat_pointinterval(
        .width = 0.95,
        mapping = aes(
            x = time,
            y = 10^log10_titer_per_ml,
            shape = detectable,
            fill = material,
            group = titer_id),
        size = 4,
        data = titer_draws,
        stroke = 2) +
    stat_pointinterval(
        .width = 0.6827,
        mapping = aes(
            x = time,
            y = 10^log10_titer_per_ml,
            shape = detectable,
            fill = material,
            group = titer_id),
        size = 8,
        fatten_point = 2,
        data = titer_draws,
        stroke = 2) +
    scale_fill_manual(values = unlist(material_colors)) +
    scale_fill_manual(values = unlist(material_colors),
                       aesthetics = "point_fill") +
    scale_color_manual(values = unlist(material_colors)) +
    shape_scale + 
    scale_y_log10_mathformat() +
    coord_cartesian(
        ylim = c(1e0, 1e6),
        xlim = c(0, 2)) +
    facet_wrap(vars(material))

# styling: no facet labels because is background plot
fit_panel <- fit_panel +
    theme_project() +
    theme(legend.position = "none") +
    xlab("time (hrs)") +
    ylab(titer_ylab)


####################################
## compose full figure from panels
####################################

cat('making full figure...\n')

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          fit_panel,
          base_width = 10,
          base_height = 8)
warnings()
