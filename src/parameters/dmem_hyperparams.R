#!/usr/bin/env Rscript

###################################
## filename: decon_hyperparams.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: contains hyperparameters for
## the Bayesian model of titer decay
## for decontamination experiments
####################################

hyperparam_list <- list(
    intercept_prior_mean = 4.5, # official dose was 10^5 TCID50
    intercept_prior_sd = 0.5,
    mode_sd_intercept = 0,
    sd_sd_intercept = 0.25,
    log_hl_prior_mean = log(0.5),
    log_hl_prior_sd = 2,
    lower_lim_decay_rate = 0,
    titer_prior_mean = 2.5,
    titer_prior_sd = 3,
    debug = FALSE)

## model options
niter <- 2000
adapt_d <- 0.95
max_tree <- 10
fixed_seed <- 888345
inits_seed <- 3462
