###-------------------------------###
### Master Thesis - Lynn Henrotte ###
###-------------------------------#############
### Finding all CYP2C9 diplotypes in GeT-RM ###
###-----------------------------------------###

## Packages ##

library("tidyverse")
library("readxl")

## Preliminaries ##

path_to_folder = paste0("C:/Users/lynnh/OneDrive/Bureaublad",
                        "/2nd Master Stat/Master Thesis/WP5/Databases")
setwd(path_to_folder)

## Obtain data ##

getrm <- read_excel("updated_getrm_calls2.xlsx")
cyp2c9_dips <- unique(getrm$CYP2C9)

## Get diplotypes ##

cyp2c9_dips2 <- list()
for (dip in cyp2c9_dips) {
  star_alleles <- strsplit(dip, split = "/")[[1]]
  allele1 <- substring(star_alleles[1], first = 2)
  allele2 <- substring(star_alleles[2], first = 2)
  if (grepl("[0-9]" , allele1) && grepl("[0-9]" , allele2)) {
    cyp2c9_dips2[[length(cyp2c9_dips2) + 1]] <- c(as.numeric(allele1),
                                                  as.numeric(allele2))
  }
}
cyp2c9_dips2

