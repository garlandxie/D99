# RLQ Analysis - Xie et al. 2018 
# Script is used to clean data for multiple databases before conducting an RLQ analysis
# Code developed by Garland Xie

#### Install packages and libraries (if you don't have them!) ####

# Load up the easypackages 
if(!require(easypackages)){
  install.packages("easypackages")
  library(easypackages)
}

# check to see if you have these packages. If not, install them
packages("magrittr")

# Load libraries
# ade4: RLQ analysis
libraries("magrittr")

# Import files --------------------------------------------------------------------

# Use relative file paths
comm <- readRDS(here("Objects", "D99_comm_matrix.rds"))
met_250 <- readRDS(here("Objects", "D99_met_250.rds"))
met_500 <- readRDS(here("Objects", "D99_met_500.rds"))
trait <- readRDS(here("Objects", "D99_trait.rds"))


# Cleaning - Community Data Matrix -------------------------------------------------

# Check to see if every species is present in the community data matrix
colnames(comm) # right number of species

# Re-order columns based on alphabetical order
comm <- comm[, order(colnames(comm))]

# Change taxon names to abbreviated versions 
# plotting purposes
taxon_change <- function (spList) {
  f <- spList %>%
    # separate genus name and species name onto a list
    strsplit(split = "_") %>%
    # create abbreviations of each species 
    # take first four letter of genus (e.g. Anthidium -> Anth)
    # take first three letters of species (e.g. manicatum -> man)
    lapply(function(x) paste(substr(x[1], 1, 3), 
                             substr(x[2], 1, 3), 
                             sep = "_"))
  # return vector of abbreviated taxon names
  return(f)
}  

# running taxon_change()
colnames(comm) <- taxon_change(colnames(comm))
rownames(trait) <- taxon_change(rownames(trait))

# small change: change Ant_man <- Anth_man
colnames(comm)[which(colnames(comm) == "Ant_man")] <- "Anth_man"
rownames(trait)[which(rownames(trait) == "Ant_man")] <- "Anth_man"

# Check to see if taxon names match in trait and community data matrix
# mismatch in number of species 
nrow(trait) # spp: 31

# Find the missing species in the trait matrix
# colnames for comm, rownames for trait
colnames(comm)[!colnames(comm) %in% rownames(trait)]

# result: Ste_ver, Coe_alt, Coe_moe, Coe_say
# basically, all kleptoparasites
# cneck previous manuscript if we removed kleptoparasites (line 276, Nov24 version)
# looks good: remove kleptoparasites from community data matrix
comm <- comm[, colnames(comm) %in% rownames(trait)]

#### Cleaning - Trait Matrix #####

# Check to see if any there any missing values in the trait database 
all(is.na(trait)) # FALSE: looks good 

# Re-order row names of traits based on alphabetical order
trait <- trait[order(rownames(trait)),]

# change trait names for plotting purposes
colnames(trait) <- c("nest.mat", "feeding.spe", "volt", 
                     "pollen.trans", "body.len", "emer.time")

# Check structure of trait matrix
# NOTE: avoid running the code twice in the same session! (i.e. NA's in trait$fee_spe)

# result: four factor variables, two numeric variables
str(trait) 

# Change the class types of each trait
# nesting material: nominal
# trait states: (1) Leaves, (2) Mud, (3) Leaves + Mud, (4) Resin, (5) Leaves Chewed, (6) None, (7) Wood Pulp
if(is.factor(trait$nes_mat)) {
  
  levels(trait$pollen.trans) <- c("Leaves", "Mud", "Leaves + Mud", 
                                  "Resin", "Leaves Chewed", "None", "Wood Pulp")
  
} else {
  
  trait$nest.mat <- trait$nest.mat %>%
                    factor(levels = c(1, 2, 3, 4, 5, 6, 7),
                           labels = c("Leaves", "Mud", "Leaves + Mud", 
                                      "Resin", "Leaves Chewed", "None", "Wood Pulp")
                         )
}

# pollen transport: nominal
# trait states: (1) Corbica, (2) Gut, (3) Scopa
if(is.factor(trait$pollen.trans)) {
  
  levels(trait$pollen.trans) <- c("Corbica", "Gut", "Scopa")
  
} else {
  trait$pollen.trans <- trait$pollen.trans %>%
                        factor(levels = c(1, 2, 3), 
                        labels = c("Corbica", "Gut", "Scopa"))
}

# feeding specialization: nominal
if(is.factor(trait$feeding.spe)) {
  
  levels(trait$feeding.spe) <- c("Polyectic", "Oligolectic")
  
} else {
  trait$feeding.spe <- trait$feeding.spe %>%
                       factor(levels = c(1,2), labels = "polyectic", "oligolectic")
}

# voltinism: quantiative 
trait$volt <- trait$volt %>%
              as.numeric()

# emergence time: quantitative
trait$emer.time <- trait$emer.time %>%
                   as.numeric()

# female body length: quantitative
trait$body.len <- trait$body.len %>%
                  as.numeric()

#Double-check trait matrix structure
str(trait) # looks good!

#### Cleaning - environmental variables ####

# Relabel colnames; current ones are too long
colnames(l_met_250) <- c("250_grass", "250_tree", "250_urban", "250_edge")
colnames(l_met_500) <- c("500_grass", "500_tree", "500_urban", "500_edge")

# Make sure landcover class types add up to 1 for both 250m and 500m scale
# function checks to see if all sites have sum of 1 for all compositional landcover classes
# arguments: 
# (1) df = data frame containing compositional landcover data types per site
# (2) vec = a vector of specific colum names (e.g. grass, tree canoy)
check_prop <- function(df, vec) { 
              # apply over rows
              f <- apply(df[, vec],  
              # sum all class types
                         MARGIN = 1, FUN = function(x) sum(x)) 
              # do all sites have a sum of 1? threshold: 0.5
              if (round(var(f), digits = 4) == 0) { 
                print("all sites have a sum of 1 for all compositional landcover classes")
              } else {
                print("something went wrong! double check your data frame")
              }
}

# test cases: 
# results: looks good!
check_prop(l_met_250, c("250_grass", "250_tree", "250_urban"))
check_prop(l_met_500, c("500_grass", "500_tree", "500_urban"))

# Re-order both environmental data matrices based on rows 
# Makes it easy to do value-matching 
l_met_250 <- l_met_250[order(rownames(l_met_250)), ]
l_met_500 <- l_met_500[order(rownames(l_met_500)), ]

# check to see if sites match between both 250m and 500m datasets
# some site labels don't match up..
rownames(l_met_250) %in% rownames(l_met_500) 

# Check number of sites per environmental data matrices 
# 250m: 167 sites, 500m: 159
# result: 8 sites missing from 500m, could be due to edge effects..
nrow(l_met_250) == nrow(l_met_500) 

# which sites? "ABB" "ADB" "ADN" "ADW" "AEN" "AEO" "AEP" "AEV"
rownames(l_met_250)[!rownames(l_met_250) %in% rownames(l_met_500)]

# discussed with N. Sookhan: perform separate analysis for each spatial scale

#### Saving data ##### 
saveRDS(comm, here("Objects", "D99_comm_matrix.rds"))
saveRDS(l_met_250, here("Objects", "D99_met_250.rds"))
saveRDS(l_met_500, here("Objects", "D99_met_500.rds"))
saveRDS(trait, here("Objects", "D99_trait.rds"))
