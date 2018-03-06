# RLQ Analysis - Xie et al. 2018 
# Script is used to create factor loadings of traits from the RLQ analysis
# Code developed by Garland Xie

# Install packages and libraries --------------------------------------------------

# Load up the easypackages 
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# check to see if you have these packages. If not, install them
packages("here")

# Load libraries
# here: constructing file paths 
libraries("here")

# Import ------------------------------------------------------------------------------

# use relative file paths
# use readRDS to prevent overwriting files
rlq_250 <- readRDS(here("Objects", "D99_rlq_250.rds"))
rlq_500 <- readRDS(here("Objects", "D99_rlq_500.rds"))

# Calculate Factor loadings -----------------------------------------------------------

# function: manually calculate factor loadings for RLQ
# arguments: (1) a "rlq dudi" R class object (so either rlq_250.rds or rlq_500.rds)
factor_loads <- function(rlq_obj) {
  
  # Manually create trait loadings for each environmental variable
  Z <- crossprod(as.matrix(rlq_obj$tab), as.matrix(rlq_obj$tab)) #cross-product matrix
  
  # Extract eigenvalues and eigenvectors from Z cross-product matrix
  eigVals <- eigen(Z)$values # eigenvalues from Z
  eigVecs <- eigen(Z)$vectors # eigenvectors from Z
  
  # create dataframe of trait loadings
  traitLoad <- data.frame(eigVecs)
  rownames(traitLoad) <- colnames(rlq_obj$tab)
  colnames(traitLoad) <- paste('Axis', 1:length(eigVals))
  
  # return trait loadings
  return(traitLoad)
  
}

# run test cases 
load_250 <- factor_loads(rlq_250)
load_500 <- factor_loads(rlq_500)

# Save files as R objects -------------------------------------------------------------

# Select first two RLQ axes 
saveRDS(load_250[, c("Axis 1", "Axis 2")], here("Objects", "D99_load_250.rds"))
saveRDS(load_500[, c("Axis 1", "Axis 2")], here("Objects", "D99_load_500.rds"))

# Write csv files ----------------------------------------------------------------------

# Select first two RLQ axess
write.csv(load_250[, c("Axis 1", "Axis 2")], here("Outputs", "D99_load_250.csv"))
write.csv(load_500[, c("Axis 1", "Axis 2")], here("Outputs" , "D99_load_500.csv"))
