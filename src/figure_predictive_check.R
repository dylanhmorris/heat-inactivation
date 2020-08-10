#!/usr/bin/env Rscript

########################################
## filename: figure_th_predictive_check.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot predictive checks for temp/humidity
## model
#######################################

suppressPackageStartupMessages(library(virusenv))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(extrafont))


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly=TRUE)
data_path <- args[1]
titer_est_data_path <- args[2]
check_results_path <- args[3]
titers_path <- args[4]
outpath <- args[5]

## read data / style files

cat("reading data (this may take a while)...\n")
check_chains <- readRDS(check_results_path)
titer_ests_chains <- readRDS(titers_path)

dat <- read_csv(data_path,
                col_types = cols())
ests_dat <- read_csv(titer_est_data_path,
                     col_types = cols())
cat("data read succesfully!\n")

## add columns where needed 
if(!"temperature" %in% names(dat)){
    if("temperature_measured" %in% names(dat)){
        dat <- dat %>%
            mutate(temperature = target_temperature)
    }
}

if(!"humidity" %in% names(dat)){
    if("target_humidity" %in% names(dat)){
        dat <- dat %>%
            mutate(humidity = target_humidity)
    }
}

if(!"drying_time" %in% names(dat)){
    dat <- dat %>%
        mutate(drying_time = 0)
}


if(!"temperature" %in% names(ests_dat)){
    if("target_temperature" %in% names(dat)){
        ests_dat <- ests_dat %>%
            mutate(temperature = target_temperature)
    }
}

if(!"humidity" %in% names(ests_dat)){
    if("target_humidity" %in% names(dat)){
        ests_dat <- ests_dat %>%
            mutate(humidity = target_humidity)
    }
}

if(!"drying_time" %in% names(ests_dat)){
   ests_dat <- ests_dat %>%
        mutate(drying_time = 0)
}

#################################
## overall plot styling
#################################
set.seed(34634) # reproducible! (since we use random draws)

gen_params <- get_params("general_param_list")
general_palette <- get_params("general_palette")
material_colors <- get_params("material_colors")

text_size = 40
detection_linesize = 0.75
titer_ylab <- expression("Virus titer (TCID"[50] * "/mL media)")
aer_ylab <- expression("Virus titer (TCID"[50] * "/L air)")
convert_mL_media_to_L_air <- 10 / 3
xlim <- c(0, max(ests_dat$time))
ylim <- c(1e-1, 1e5)
LOD_log10 <- 0.5
LOD <- 10^LOD_log10
line_alpha <- gen_params[["line_alpha"]]
interval_grey <- general_palette[["interval_grey"]]

conversion_factor <- 1
n_lines <- 10


##################################################
## calculate posterior draws for regression lines
##################################################


cat("extracting draws for decay rates / intercepts (this may also take a while)...\n")

tidy_draws <- check_chains %>%
    spread_draws(true_titers_pred[titer_id])


cat("extracting positive wells...\n")
pos_wells <- ests_dat %>%
    group_by(titer_id) %>%
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))

print(pos_wells)

tidy_draws <- dat %>%
    distinct(titer_id,
             .keep_all = TRUE) %>%
    select(trial_unique_id,
           virus,
           material,
           temperature,
           humidity,
           drying_time,
           titer_id,
           time,
           treatment) %>% 
    inner_join(tidy_draws,
               by = "titer_id") %>%
    filter(time >= drying_time) %>%
    mutate(time = time - drying_time)

cat('extracting titer estimates...\n')

titer_ests_draws <- titer_ests_chains %>%
    spread_draws(log10_tcid50[titer_id])


## get human readable names and detectability
titer_ests_draws <- ests_dat %>%
    distinct(titer_id,
             .keep_all = TRUE) %>%
    select(replicate,
           virus,
           time,
           material,
           drying_time,
           temperature,
           titer_id,
           humidity,
           treatment) %>% 
    inner_join(titer_ests_draws,
               by = "titer_id") %>%
    inner_join(pos_wells,
               by = "titer_id") %>%
    mutate(detectable = n_pos > 1) %>%
    filter(material == "Plastic")

## check that small and undetectable titers have
## higher variance if there are fewer wells

cat("converting units and setting undetectables to detection limit\n")

## convert from TCID50/(0.1mL) to TCID50/mL
## and visualize 0 positive well titers at
## the traditional LOD
titer_ests_draws <- titer_ests_draws %>%
    mutate(log10_tcid50 = ifelse(
               detectable,
               log10_tcid50 + 1,
               LOD_log10)) %>%
    arrange(desc(time), desc(material)) ## so masks are on top 

print(titer_ests_draws)

## adding one to convert to per mL from per 0.1 mL

shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))

plot_dat <- titer_ests_draws %>%
    filter(time > drying_time) %>%
    mutate(time = time - drying_time)

scale <- 2 * sqrt(xlim[2])
if(grepl("posterior", outpath))
    scale <- xlim[2]

panel <- tidy_draws %>%
    ggplot(aes(x = time,
               y = 10^(1 + true_titers_pred) * conversion_factor,
               color = material,
               fill = material)) +
    geom_hline(aes(yintercept = LOD),
               size = 2,
               linetype = "dotted") +
    stat_eye(color = "black",
             scale = scale,
             slab_color = "black",
             size = 15) +
    stat_pointinterval(
        .width = 0.95,
        mapping = aes(x = time,
                      y = 10^log10_tcid50 * conversion_factor,
                      shape = detectable,
                      fill = material,
                      group = titer_id),
        data = plot_dat,
        point_size = 6,
        size = 7,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9) +
    stat_pointinterval(
        .width = 0.6827,
        mapping = aes(x = time,
                      y = 10^log10_tcid50 * conversion_factor,
                      shape = detectable,
                      fill = material,
                      group = titer_id),
        data = plot_dat,
        point_size = 6,
        size = 14,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9) +
    scale_fill_manual(values = unlist(material_colors)) +
    scale_fill_manual(values = unlist(material_colors),
                       aesthetics = "point_fill") +
    scale_color_manual(values = unlist(material_colors)) +
    shape_scale + 
    scale_y_log10_mathformat(expand = c(0, 0)) +
    coord_cartesian(ylim = ylim,
                    xlim = xlim) +
    facet_grid(rows = vars(humidity),
               cols = vars(temperature))

## styling: no facet labels because is background plot
panel <- panel +
    theme_project(base_size = text_size) +
    xlab("Time (hrs)") +
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
