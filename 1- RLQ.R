# RLQ Analysis - Xie et al. 2018 
# Script is used to conduct an RLQ analysis to determine trait-environment relationships
# Code developed by Garland Xie

# Install packages and libraries --------------------------------------------------

# Load up the easypackages 
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# check to see if you have these packages. If not, install them
packages("ade4", "here")

# Load libraries
# ade4: RLQ analysis
libraries("ade4", "here")

# Import Files --------------------------------------------------------------------

# NOTE: be careful of overwriting existing files in global env
# Use relative file paths

comm <- readRDS(here("Objects", "D99_comm_matrix.rds"))
met_250 <- readRDS(here("Objects", "D99_met_250.rds"))
met_500 <- readRDS(here("Objects", "D99_met_500.rds"))
trait <- readRDS(here("Objects", "D99_trait.rds"))

# RLQ analysis ---------------------------------------------------------------------


# arguments:
# (1) R - site x environment
# (2) Q - species x traits
# (3) L - site x species
# assume that you have consistent sites, species and traits 
rlq_custom <- function(dudiR, dudiQ, dudiL) {
  
  # Make sure sites match between R and L
  dudiL <- dudiL[rownames(dudiL) %in% rownames(dudiR), ]
  
  # Make sure species match between Q and L 
  # excludes kleptoparasites
  dudiL <- dudiL[, colnames(dudiL) %in% rownames(dudiQ)]
  
  # making the L Tabe (correspondence analysis)
  L <- dudi.coa(dudiL, scannf = F)
  
  # making the R Table (principal component analysis), constrained row weights by species abundance
  R <- dudi.pca(dudiR, row.w = L$lw, scannf = F)
  
  # making the Q Table (hillsmith analysis), constrained column weights by species abundance
  # alternatively, dudi.mix works for ordinal variables, but no options for row weights.. 
  Q <- dudi.hillsmith(dudiQ, row.w = L$cw, scannf = F)
  
  # Combining R, Q and Ltables
  RLQ <- rlq(dudiR = R, dudiQ = Q, dudiL = L)
  
  # return RLQ 
  RLQ
}

# RLQ analysis for 250m and 500m scale
# Choose the first two axes since they explain the highest amount of variation
rlq_250 <- rlq_custom(dudiR = met_250, dudiQ = trait, dudiL = comm)
rlq_500 <- rlq_custom(dudiR = met_500, dudiQ = trait, dudiL = comm)


# Saving files --------------------------------------------------------------

saveRDS(rlq_250, here("Objects", "D99_rlq_250.rds"))
saveRDS(rlq_500, here("Objects", "D99_rlq_500.rds"))
