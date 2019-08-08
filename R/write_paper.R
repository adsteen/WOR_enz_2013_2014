# Rewrite of 2016_06_03_full_WOR_peptidase

library(tidyverse)
library(broom)
library(lmstats)
library(scales)
library(cowplot)

source("R/predict_mm.R") # May not be necessary if I use purrr properly
source("R/safe_coef.R")
source("R/mutate_cond.R")
source("R/data_lm.R")
source("R/get_slope.R")

print.plots <- TRUE
print.deep.cuts <- FALSE
save.plots <- FALSE

# Set a visual theme for all the plots
theme_set(theme_bw() + 
            theme(panel.grid.major=element_blank(), 
                  panel.grid.minor=element_blank(),
                  text=element_text(size=7.5),
                 plot.margin=grid::unit(c(0,0,0,0), "mm")))

# data for approx depth of smtz
smtz <- data.frame(quantDepth=c(40, 60), ymin=c(0, 0), ymax=c(1, 1), rel.Vmax=c(1, 1))

# Correspondence between categorical depth (1:5) and quantitative depth (cmbsf)
quant_depth <- data.frame(depth = 1:5, quant.depth = c(1.5, 4.5, 28.5, 58.5, 82.5))

######
# Read in raw data, calculate kinetics and make Fig 1
######
source("R/2018_ms_read_and_fit_data.R")

####
# Make Km vs depth figure
####

source("2018_p_Km_v_depth.R")


# Make the plot of respiration rate per cell
source("R/2018_vmax_rel_G.R")

# Gene stuff
source("R/2019_07_10_phyla_abundances.R")

source("R/2018_11_28_Baker_reanalysis.R")

# print(p_Vmax.cell.norm)

# Calculate Vmax relative to G

