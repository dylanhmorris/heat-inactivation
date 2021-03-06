#!/usr/bin/env Rscript

###################################
## shared styling for all plots
##
####################################

## styling for text

power10_breaks <- scales::trans_breaks("log10",
                                       function(x) 10^x)
power10_format <- scales::trans_format("log10",
                                       scales::math_format(10^.x))


## virus order

#' virus_order
#'
#' order in which to label viruses
#'
virus_order <- c("SARS-CoV-2",
                 "SARS-CoV-1",
                 "MERS-CoV",
                 "HCoV-229E",
                 "HCoV")

#' scale_x_log10_mathformat
#' 
#' convenience scale function
#' for putting x breaks on
#' powers of 10 in ggplot
#'
#' @param ... arguments passed to scale_x_continuous()
#'
#' @return a scale function
#' 
#'@export
scale_x_log10_mathformat <- function(...){
    return (ggplot2::scale_x_continuous(trans = "log10",
                                        breaks = power10_breaks,
                                        label = power10_format,
                                        ...))
}

#'
#'
#' convenience scale function
#' for putting y breaks on
#' powers of 10 in ggplot
#'
#' @param ... arguments passed to scale_y_continuous()
#'
#' @return a scale function
#' 
#'@export
scale_y_log10_mathformat <- function(...){
    return (ggplot2::scale_y_continuous(trans = "log10",
                                        breaks = power10_breaks,
                                        label = power10_format,
                                        ...))
}

param_list <- list(

    ## display names for viruses
    "virus_display_names" = list(
        "sars2" = "SARS-CoV-2",
        "sars" = "SARS-CoV-1",
        "mers" = "MERS-CoV",
        "flu" = "Influenza"),

    ## color palette
    "general_palette" = list(
        "interval_grey" = "#525252"
    ),
    
    "virus_colors" = list(
        "SARS-CoV-2" = "#7c0885",
        "SARS-CoV-1" = "#407dc2",
        "MERS-CoV" = "#9bc2b2",
        "HCoV" = "#e6e4c5",
        "HCoV-229E" = "#e6e4c5",
        "MHV" = "#518591",
        "TGEV" = "#7570b3"),

    "material_colors" = list(
        "N95 mask" = "#554fa8",
        "Plastic" = "#7c0885",
        "Aerosols" = "#54e3ff",
        "Cardboard" = "#615242",
        "Copper" = "#da8a67",
        "Steel" = "#adbab8",
        "DMEM" = "#609c70",
        "Nasal_wet" = "#609c70",
        "Nasal_dry" = "#609c70",
        "Sputum_wet" = "#609c70",
        "Sputum_dry" = "#609c70",
        "Uncovered plate oven" = "#609c70",
        "Closed vial heat block" = "#609c70",
        "Closed vial oven" =  "#609c70",
        "Covered plate oven" =  "#609c70"),

    "fluid_type_colors" = list(
        "Nasal" = "#609c70",
        "Sputum" = "#264475"),

    "treatment_colors" = list(
        "Control" = "grey",
        "Ethanol" = "#d2db1f",
        "Heat" = "#ff8b0f",
        "VHP" = "#0fc7ff",
        "UV" = "#a503ab"),

    "model_colors" = list(
        "Mechanistic" = "green",
        "Null" = "gray"),

    "general_param_list" = list(
        "panel_fig_height" = 20,
        "panel_fig_aspect" = 2.5,
        "line_fineness" = 500,  ## fineness of line plotting
        "n_lines" = 50,    ## number of random posterior lines to plot
        "line_alpha" = 0.1,  ## transparency of overplotted lines
        "line_size" = 2,
        "detection_linestyle" = "dashed")
)


## styling for figures

#' get_params
#'
#' get the desired list of plotting
#' parameters
#' 
#' @param parameter_list_name name of the parameter list to get, one of: \cr
#'
#' virus_display_names: list of display names for studied viruses \cr
#'
#' general_palette: general attributes of the project color palette \cr
#'
#' virus_colors: list of colors associated to each studied virus \cr
#'
#' material_colors: list of colors associated to each material
#' on which decay was studied \cr
#'
#' fluid_type_colors: list of colors associated with human fluids
#' of interest \cr
#'
#' treatment_colors: list of colors associated with various
#' decontamination treatments \cr
#'
#' general_param_list: general plot styling parameters, such
#' as line weight and style \cr 
#' 
#' @export
get_params <- function(parameter_list_name){


    return( param_list[[parameter_list_name]] )
}


#' read_data_for_plotting
#'
#' function to read in cleaned data and format
#' it as a tibble that can be used for plotting
#'
#' @param data_path path to data to read in
#' @return tibble of formatted data
#' @export
read_data_for_plotting <- function(data_path) {
    dat <- readr::read_csv(data_path,
                           col_types = cols())

    ## add unlogged titers for plotting
    if ("log10_titer" %in% names(dat) ) {
        dat$titer = 10^dat$log10_titer
    }
   
    ## make sure the viruses/strains display
    ## in our order or interest (set below)
    dat$virus <- factor(dat$virus,
                        virus_order)
   
    if ("strain" %in% names(dat) ) {
        dat$strain <- factor(dat$strain,
                             strain_order)
    }
    
    return (dat)
}


vector_formats <- c(".pdf", ".svg")

#' save_plot
#'
#' override cowplot's save_plot
#' to include font embedding
#'
#' @param outpath path to save the figure to
#' @param fig ggplot figure to save
#' @param ... other parameters passed to cowplot::save_plot()
#' @return nothing
#' @export
save_plot <- function(outpath, fig, ...){
    cowplot::save_plot(outpath,
                       fig,
                       ...)
    if(any(endsWith(outpath, vector_formats)))
        extrafont::embed_fonts(file = outpath)
}

#' theme_project
#'
#' variant of theme_classic() ggplot theme for
#' this project
#' 
#' @param base_size base font size for the theme, passed to
#' theme_classic (default 30)
#' @param ... other parameters passed to theme_classic()
#' @return theme
#' @export
theme_project <- function(base_size = 30,
                          ...){
    x_axis_margin <- ggplot2::margin(t = base_size / 2)
    y_axis_margin <- ggplot2::margin(r = base_size / 2)

    ggplot2::theme_classic(
                 base_size = base_size,
                 ...) +
    cowplot::background_grid(major = "xy",
                             minor = "none",
                             size.major = 0.5) +
    ggplot2::theme(axis.title.x = element_text(
                       margin = x_axis_margin),
                   axis.title.y = element_text(
                       margin = y_axis_margin),
                   strip.background = element_blank(),
                   panel.border = element_rect(color = "black",
                                               fill = NA,
                                               size = 2),
                   text = element_text(family = "CM Sans"))

}


#' theme_default_legend
#'
#' default legend styling for this project
#' 
#' @return legend theme
#' @export
theme_default_legend <- function(){

    ggplot2::theme(legend.margin = margin(t = -10,
                                          l = 0,
                                          r = 5,
                                          b = 0),
                   legend.title = element_blank(),
                   legend.box.background = element_rect(color = "black",
                                                        size = 2))
}
