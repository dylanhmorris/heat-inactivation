#####################################
# name: Makefile
# author: Dylan Morris <dhmorris@princeton.edu>
#
# Makefile to generate analyses
# for Gamble et al study of virus
# heat inactivation as a function
# of heating procedure
####################################

#####################################
# Directory structure
####################################

default: all

SRC := src
OUTPUT := out
DATA := dat
RAW := $(DATA)/raw
CLEANED := $(DATA)/cleaned
PARAMS := $(SRC)/parameters

MCMC_CHAINS := $(OUTPUT)/mcmc_chains
FIGURES = $(OUTPUT)/figures
TABLES = $(OUTPUT)/tables

#####################################
# File extensions and the like
####################################
CHAINS_SUFFIX = _chains.Rds
PRIOR_CHECK_SUFFIX = _prior_check.Rds

#####################################
# Expected bash settings
#
# Check these vs your local
# machine setup if you are having
# difficulty reproducing the
# analysis
#####################################

MKDIR := @mkdir -p
RM := rm -rf

R_OPTIONS = --vanilla
R_COMMAND := Rscript $(R_OPTIONS)

# R creates a blank Rplots.pdf when run
# from the command line to produce
# a figure. This removes it.
FIG_CLEANUP = @$(RM) Rplots.pdf

#####################################
# Installation / dependencies
#
# Rules for prepping analysis
#####################################

.PHONY: depend

depend:
	$(R_COMMAND) $(SRC)/install_needed_packages.R


################################
# model names
################################

DMEM_HL_NAME = dmem-hl
DMEM_TITER_NAME = dmem-titer
MODEL_NAMES = $(DMEM_HL_NAME) $(DMEM_TITER_NAME)

#####################################
# data locations
#####################################
DMEM_DATAFILE = dmem_data.csv

DATAFILES = $(DMEM_DATAFILE)

CLEANED_DATA = $(addprefix $(CLEANED)/, $(DATAFILES))

###########################
## associate data to models

$(DMEM_HL_NAME)_FITTING_DATA = $(CLEANED)/$(DMEM_DATAFILE)
$(DMEM_TITER_NAME)_FITTING_DATA = $(CLEANED)/$(DMEM_DATAFILE)

#####################################
# code locations
#####################################

## scripts for cleaning and fitting
DEFAULT_CLEANING_SCRIPT = $(SRC)/clean_data.R
DEFAULT_FITTING_SCRIPT = $(SRC)/fit_model.R
DIAGNOSTIC_SCRIPT = $(SRC)/chain_diagnostics.R

## associate scripts to model names
$(DMEM_HL_NAME)_FITTING_SCRIPT = $(DEFAULT_FITTING_SCRIPT)
$(DMEM_TITER_NAME)_FITTING_SCRIPT = $(DEFAULT_FITTING_SCRIPT)

## stan model files
DEFAULT_HL_SRC = $(SRC)/well_decay_model.stan
DEFAULT_INFER_SRC = $(SRC)/well_titer_estimates.stan

## associate stan files to model names
$(DMEM_HL_NAME)_MODEL_SRC = $(DEFAULT_HL_SRC)
$(DMEM_TITER_NAME)_MODEL_SRC = $(DEFAULT_INFER_SRC)

#####################################
# parameter locations
#####################################

## paths to hyperparameters
DEFAULT_HYPERS = $(PARAMS)/dmem_hyperparams.R

## associate params to model names
$(DMEM_HL_NAME)_HYPERS = $(DEFAULT_HYPERS)
$(DMEM_TITER_NAME)_HYPERS = $(DEFAULT_HYPERS)

#####################################
# mcmc output locations
#####################################

## fit chains
FIT_CHAINS = $(addprefix $(MCMC_CHAINS)/, $(addsuffix $(CHAINS_SUFFIX), $(MODEL_NAMES)))

## prior predictive checks

PRIOR_CHECK_NAMES = $(addsuffix $(PRIOR_CHECK_SUFFIX), $(MODEL_NAMES))
PRIOR_CHECK_CHAINS = $(addprefix $(MCMC_CHAINS)/, $(PRIOR_CHECK_NAMES))

CHAIN_PATHS = $(PRIOR_CHECK_CHAINS) $(FIT_CHAINS)

CHAIN_DIAGNOSTICS = $(OUTPUT)/chain_diagnostics.csv

#####################################
# Rules
#
# definition of dependency
# tree and specification of
# rules for doing stuff
#####################################

##########################
# rules for data cleaning
##########################
$(CLEANED)/%_data.csv: $(DEFAULT_CLEANING_SCRIPT) $(RAW)/%_data.csv

	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $^ $@

#####################################
# rules for model fitting and post-processing
#####################################
.SECONDEXPANSION:

# generic model fitting rule
$(MCMC_CHAINS)/%$(CHAINS_SUFFIX): $$($$*_FITTING_SCRIPT) $$($$*_MODEL_SRC) $$($$*_FITTING_DATA) $$($$*_HYPERS)
	@echo Making $@... "\n\n"
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@ FALSE TRUE

$(MCMC_CHAINS)/%$(PRIOR_CHECK_SUFFIX): $$($$*_FITTING_SCRIPT) $$($$*_MODEL_SRC) $$($$*_FITTING_DATA) $$($$*_HYPERS)
	@echo "\nMaking" $@... "\n"
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@ TRUE TRUE

$(CHAIN_DIAGNOSTICS): $(DIAGNOSTIC_SCRIPT) $(CHAIN_PATHS)
	$(MKDIR) $(OUTPUT)
	$(R_COMMAND) $^ $@

#####################################
# rules for figures
#####################################

$(FIGURES)/%.pdf: $(SRC)/%.R $(CLEANED)/$(DMEM_DATAFILE) $(FIT_CHAINS)
	$(MKDIR) $(FIGURES)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


## predictive checks

$(FIGURES)/figure-%-prior-check.pdf: $(SRC)/figure-predictive-check.R $$($$*_FITTING_DATA) $(MCMC_CHAINS)/$$*$(PRIOR_CHECK_SUFFIX) $(MCMC_CHAINS)/$(DMEM_TITER_NAME)$(CHAINS_SUFFIX)
	$(MKDIR) $(FIGURES)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURES)/figure-%-posterior-check.pdf: $(SRC)/figure-predictive-check.R $$($$*_FITTING_DATA) $(MCMC_CHAINS)/$$*$(CHAINS_SUFFIX) $(MCMC_CHAINS)/$(DMEM_TITER_NAME)$(CHAINS_SUFFIX)
	$(MKDIR) $(FIGURES)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)


PRIOR_CHECK_FIGS = $(addprefix figure-, $(addsuffix -prior-check.pdf, $(MODEL_NAMES)))

POSTERIOR_CHECK_FIGS = $(addprefix figure-, $(addsuffix -posterior-check.pdf, $(MODEL_NAMES)))

FIGURE_NAMES = figure-dmem-decay.pdf figure-dmem-halflives.pdf $(PRIOR_CHECK_FIGS) $(POSTERIOR_CHECK_FIGS)

FIGURE_PATHS = $(addprefix $(FIGURES)/, $(FIGURE_NAMES))

#####################################
# rules for tables
#####################################

$(TABLES)/%.tex: $(SRC)/%.R $(CLEANED)/$(DMEM_DATAFILE) $(FIT_CHAINS)
	$(MKDIR) $(TABLES)
	$(R_COMMAND) $^ $@

TABLE_NAMES = table-dmem-halflives.tex
TABLE_PATHS = $(addprefix $(TABLES)/, $(TABLE_NAMES))

#####################################
# convenience rules for making
# various quantities
#####################################
.PHONY: data
data: $(CLEANED_DATA)

.PHONY: chains
chains: $(CHAIN_PATHS)

.PHONY: diagnostics
diagnostics: $(CHAIN_DIAGNOSTICS)

.PHONY: figures
figures: $(FIGURE_PATHS) 

.PHONY: tables
tables: $(TABLE_PATHS)

.PHONY: echo_figures echo_chains

echo_figures:
	@echo
	@echo $(FIGURE_PATHS)
	@echo

echo_chains:
	@echo
	@echo $(CHAIN_PATHS)
	@echo

## remove emacs tempfiles, etc.
.PHONY deltemp:
	$(RM) *~*
	$(RM) $(SRC)/*~*
	$(RM) $(SRC)/*#*
	$(RM) $(PARAMS)/*~*
	$(RM) $(PARAMS)/*#*

.PHONY: clean
clean: deltemp
	$(RM) $(OUTPUT)
	$(RM) $(CLEANED)

all: depend data chains diagnostics figures tables
