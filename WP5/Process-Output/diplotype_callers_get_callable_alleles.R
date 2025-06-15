### Master Thesis: synthetic read simulation for PGx ###
### Lynn Henrotte 2024-2025 ###

### Packages ###

library("tidyverse")
library("ggplot2")
library("ggpubr")

###---------------------------###
###-- Find callable alleles --###
###---------------------------###

# Read in data
callable_alleles <- read.csv("Callable_alleles.csv", sep = "\t")

# Get all genes to be included
genes <- c("CYP2C9", "CYP2C19")
callable_alleles_sub <- callable_alleles[, colnames(callable_alleles) %in% 
                                           genes]

genes <- colnames(callable_alleles)[2:15]

# Get caller names 
callers <- callable_alleles$Caller

# Initiate lists
callable_alleles_list <- list()

# Gather callable star alleles in a list of lists (one per gene)
for (gene in genes) {
  
  # Define list to store callable alleles for the gene as vectors
  temp_list <- list()
  
  for (row in 1:6) {
    
    # Obtain data
    data_row <- callable_alleles[gene][row,]
    data_row <- substr(data_row, 3, nchar(data_row) - 2) # Exclude [' and '] parts
    
    # Split data to obtain vector of star alleles
    alleles <- stringr::str_split(data_row, pattern = "', '")[[1]]
    
    # If no star alleles can be called, set list entry to NA
    if (alleles[1] == "") {
      alleles <- NA
    }
    
    # Add vector (or NA) to list
    temp_list[[row]] <- alleles
  }
  
  # Set list item names to names of callers
  names(temp_list) <- callers
  
  # Remove NA entries
  temp_list2 <- temp_list[!is.na(temp_list)]
  
  # Add list for gene to overall list
  callable_alleles_list[[ gene ]] <- temp_list2
  
}

all_cyp2c9_alleles <- 1:85
all_cyp2c19_alleles <- c(1:19,22:26,28:39)

#### PyPGx #####

## CYP2C9 ##

# Extract information on CYP2C9 and PyPGx
pypgx_cyp2c9_callable_alleles <- callable_alleles_list[["CYP2C9"]][["PyPGx"]]

# Convert star alleles to numbers for easier representation + sorting
pypgx_cyp2c9_callable_num <- as.numeric(substring(pypgx_cyp2c9_callable_alleles, 2, 3))

# Sort the callable CYP2C9 alleles
pypgx_cyp2c9_callable_num <- sort(pypgx_cyp2c9_callable_num)

# Get non-callable CYP2C9 alleles
setdiff(all_cyp2c9_alleles, pypgx_cyp2c9_callable_num)

## CYP2C19 ##

# Extract information on CYP2C19 and PyPGx
pypgx_cyp2c19_callable_alleles <- callable_alleles_list[["CYP2C19"]][["PyPGx"]]

# Convert star alleles to numbers for easier representation + sorting
pypgx_cyp2c19_callable_num <- as.numeric(substring(pypgx_cyp2c19_callable_alleles, 2, 3))

# Sort the callable CYP2C19 alleles
pypgx_cyp2c19_callable_num <- sort(pypgx_cyp2c19_callable_num)

# Get non-callable CYP2C19 alleles
setdiff(all_cyp2c19_alleles, pypgx_cyp2c19_callable_num)

#### Aldy #####

## CYP2C9 ##

# Extract information on CYP2C9 and Aldy
aldy_cyp2c9_callable_alleles <- callable_alleles_list[["CYP2C9"]][["Aldy"]]

# Convert star alleles to numbers for easier representation + sorting
aldy_cyp2c9_callable_num <- as.numeric(substring(aldy_cyp2c9_callable_alleles, 2, 3))

# Sort the callable CYP2C9 alleles
aldy_cyp2c9_callable_num <- sort(aldy_cyp2c9_callable_num)

# Get non-callable CYP2C9 alleles
setdiff(all_cyp2c9_alleles, aldy_cyp2c9_callable_num)

## CYP2C19 ##

# Extract information on CYP2C19 and Aldy
aldy_cyp2c19_callable_alleles <- callable_alleles_list[["CYP2C19"]][["Aldy"]]

# Convert star alleles to numbers for easier representation + sorting
aldy_cyp2c19_callable_num <- as.numeric(substring(aldy_cyp2c19_callable_alleles, 2, 3))

# Sort the callable CYP2C19 alleles
aldy_cyp2c19_callable_num <- sort(aldy_cyp2c19_callable_num)

# Get non-callable CYP2C19 alleles
setdiff(all_cyp2c19_alleles, aldy_cyp2c19_callable_num)

#### StellarPGx #####

## CYP2C9 ##

# Extract information on CYP2C9 and StellarPGx
stellarpgx_cyp2c9_callable_alleles <- callable_alleles_list[["CYP2C9"]][["StellarPGx"]]

# Convert star alleles to numbers for easier representation + sorting
stellarpgx_cyp2c9_callable_num <- as.numeric(substring(stellarpgx_cyp2c9_callable_alleles, 2, 3))

# Sort the callable CYP2C9 alleles
stellarpgx_cyp2c9_callable_num <- sort(stellarpgx_cyp2c9_callable_num)

# Get non-callable CYP2C9 alleles
setdiff(all_cyp2c9_alleles, stellarpgx_cyp2c9_callable_num)

## CYP2C19 ##

# Extract information on CYP2C19 and StellarPGx
stellarpgx_cyp2c19_callable_alleles <- callable_alleles_list[["CYP2C19"]][["StellarPGx"]]

# Convert star alleles to numbers for easier representation + sorting
stellarpgx_cyp2c19_callable_num <- as.numeric(substring(stellarpgx_cyp2c19_callable_alleles, 2, 3))

# Sort the callable CYP2C19 alleles
stellarpgx_cyp2c19_callable_num <- sort(stellarpgx_cyp2c19_callable_num)

# Get non-callable CYP2C19 alleles
setdiff(all_cyp2c19_alleles, stellarpgx_cyp2c19_callable_num)
