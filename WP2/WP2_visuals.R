#####################
### MASTER THESIS ###
#####################

### WP2: assessment of the extent of missing PGx variation ###
### Lynn Henrotte

##----------------------------------------------------------------------------##

### PACKAGES ###

library("ggplot2")
library("tidyverse")
library("ggpubr")
library("ComplexUpset")
library("VennDiagram")
library("ursaPGx")
library("stringr")

##----------------------------------------------------------------------------##

# Set working directory
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP2")

##----------------------------------------------------------------------------##

# Collect all star allele definitions on PharmVar (version 6.2.3)

defined_alleles <- read.csv("Defined_alleles_all.csv", sep = "\t")
defined_alleles <- defined_alleles %>% arrange(Gene)
genes <- defined_alleles$Gene

# Initiate list
defined_alleles_list <- list()

# Gather defined star alleles in a list of lists (one per gene)
for (gene in genes) {
  
  # Extract all PharmVar recognized star alleles for gene
  data_all <- defined_alleles[defined_alleles$Gene == gene, ]$Alleles
  
  # Exclude [' and '] parts
  data_all <- substr(data_all, 3, nchar(data_all) - 2) 
  
  # Split string to obtain individual alleles
  alleles_all <- stringr::str_split(data_all, pattern = "', '")[[1]]
  
  # Add to list
  defined_alleles_list[[ gene ]] <- alleles_all
  
}

##----------------------------------------------------------------------------##

# Load in allele coverage results
cov_callers <- read.csv("Diplotype_Caller_Allele_Coverage.csv", sep = "\t")
cov_callers <- cov_callers %>% arrange(Gene)
cov_callers$ursaPGx_Coverage <- ursaPGx_coverages

cov_databases <- read.csv("Database_Allele_Coverage.csv", sep = "\t")
cov_databases <- cov_databases %>% arrange(Gene)

# Move gene names to row names
row.names(cov_callers) <- cov_callers$Gene
row.names(cov_databases) <- cov_databases$Gene
cov_callers <- cov_callers[, 2:7]
cov_databases <- cov_databases[, 2:3]

# Compute mean among non-zero coverages:
mean_cov_callers <- apply(X = cov_callers, MARGIN = 2, FUN = mean, na.rm = TRUE)

# Replace NA's by 0% coverage
cov_callers[is.na(cov_callers)] <- 0
cov_databases[is.na(cov_databases)] <- 0

## Compute average coverage of callers including zeroes
mean_cov_callers2 <- apply(X = cov_callers, MARGIN = 2, FUN = mean)
apply(X = cov_callers, MARGIN = 2, FUN = median)

##################################
### DATABASES: ALLELE COVERAGE ###
##################################

## Prepare data

gene_names_db <- row.names(cov_databases)
genes_db <- rep(gene_names_db, times = 2)

database <- rep(c("GeT-RM", "Star Allele Search"), each = length(gene_names_db))

coverage_db <- c(cov_databases$GeT.RM, cov_databases$Star.allele.search)

labs_db <- scales::percent(coverage_db, accuracy = 1)

data_db <- data.frame("Gene" = genes_db, "Database" = database, "Coverage" = coverage_db,
                      "Labels"= labs_db)

## Create plot

p_db <- ggplot(data = data_db, aes(fill = Database, color = Database, y = Coverage, 
                                   x = Gene, label = Labels)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) + theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5,
            fontface = "bold", show.legend = FALSE) + labs(y = "Allele Coverage") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold.italic"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 18, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 20))  +
  ggtitle("Star allele coverage of databases for different genes") +
  scale_fill_manual(values = c("royalblue2", "maroon3")) +
  scale_color_manual(values = c("royalblue4", "maroon"))

## Create plot for CYP2C9 and CYP2C19 only

data_db_sub <- data_db[data_db$Gene %in% c("CYP2C9", "CYP2C19"),]

p_db_sub <- ggplot(data = data_db_sub, 
                   aes(fill = Database, color = Database, y = Coverage, 
                       x = Gene, label = Labels)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) + theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 5.5,
            fontface = "bold", show.legend = FALSE) + labs(y = "Allele Coverage") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold.italic"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 16, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 20, hjust = 0.5))  +
  guides(fill = guide_legend(nrow=2)) +
  ggtitle("Star allele coverages of PGx databases\nfor CYP2C9 and CYP2C19") +
  scale_fill_manual(values = c("royalblue2", "maroon3")) +
  scale_color_manual(values = c("royalblue4", "maroon"))

############################################
### STAR ALLELE CALLERS: ALLELE COVERAGE ###
############################################

## ursaPGx: find star allele coverages

genes <- defined_alleles$Gene

ursaPGx_callable_alleles <- ursaPGx::availableHaplotypes(build = "GRCh38")
ursaPGx_callable_alleles_list <- list()
ursaPGx_coverages <- c()

for (gene in genes) {
  
  # ursaPGx callable alleles for the gene
  alleles_ursaPGx <- ursaPGx_callable_alleles[grep(gene, ursaPGx_callable_alleles)]
  alleles_ursaPGx <- str_extract(alleles_ursaPGx, "\\*\\d+")
  
  # all defined alleles for the gene
  alleles_defined <- defined_alleles_list[[ gene ]]
  
  # only retain ursaPGx callable alleles that are defined on PharmVar
  alleles_ursaPGx <- intersect(alleles_ursaPGx, alleles_defined)
  
  # Compute coverage (+1 because the wild type is never listed)
  ursaPGx_coverages <- c(ursaPGx_coverages, 
                         (length(alleles_ursaPGx) + 1)/length(alleles_defined))
  
  # Save ursaPGx callable alleles for the gene
  ursaPGx_callable_alleles_list[[ gene ]] <- alleles_ursaPGx
  
}

##----------------------------------------------------------------------------##

## Prepare data

gene_names_call <- row.names(cov_callers)
genes_call <- rep(gene_names_call, times = 6)

callers <- rep(c("PyPGx", "StellarPGx", "PharmCAT", "ursaPGx", "PAnno", "Aldy"), 
               each = length(gene_names_call))

coverage_call <- c(cov_callers$PyPGx_Coverage, cov_callers$StellarPGx_Coverage,
              cov_callers$PharmCAT_Coverage, cov_callers$ursaPGx_Coverage,
              cov_callers$PAnno_Coverage, cov_callers$Aldy_Coverage)

labs_call <- scales::percent(coverage_call, accuracy = 1)

data_call <- data.frame("Gene" = genes_call, "Caller" = callers, 
                        "Coverage" = coverage_call, "Labels"= labs_call)

## Create plot

p_call <- ggplot(data = data_call, aes(fill = Caller, y = Coverage, x = Gene,
                           label = Labels)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") + theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(y = "Allele Coverage", fill = "Star Allele Caller") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold.italic"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 16, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 20)) +
  guides(fill = guide_legend(nrow = 1)) +
  ggtitle("Star allele coverages of six diplotype callers for different genes")

ggarrange(p_db, p_call, ncol = 1, nrow = 2)

## Create plot for CYP2C9 and CYP2C19 specifically

data_call_sub <- data_call[data_call$Gene %in% c("CYP2C9", "CYP2C19"),]

p_call_sub <- ggplot(data = data_call_sub, 
                     aes(fill = Caller, y = Coverage, x = Gene,
                         label = Labels)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.8, color = "black") + theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(y = "Allele Coverage", fill = "Star Allele Caller") +
  geom_text(aes(color = Caller), position = position_dodge(width = 0.8), vjust = -0.5, size = 5,
            fontface = "bold", show.legend = FALSE) +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold.italic"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 16, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 20, hjust = 0.5)) +
  guides(fill = guide_legend(nrow = 2)) +
  ggtitle("Star allele coverages of 5 diplotype callers\nfor CYP2C9 and CYP2C19")

ggarrange(p_db_sub, p_call_sub, ncol = 2)

###############################
### DATABASES: VENN DIAGRAM ###
###############################

# Read in data
covered_alleles <- read.csv("Covered_alleles_databases.csv", sep = "\t")
covered_alleles$Database <- c("GeT-RM", "StarAlleleSearch")

genes <- colnames(covered_alleles)[-c(1)]

databases <- covered_alleles$Database

# Initiate lists
covered_alleles_list <- list()

# Gather callable star alleles in a list of lists (one per gene)
for (gene in genes) {
  
  # Define list to store covered alleles for the gene as vectors
  temp_list <- list()
  
  for (row in 1:2) {
    
    # Obtain data
    data_row <- covered_alleles[gene][row,]
    
    if (data_row == "[]") {
      temp_list[[row]] <- NA
    }
    else {
      data_row <- substr(data_row, 3, nchar(data_row) - 2) # Exclude [' and '] parts
      
      # Split data to obtain vector of star alleles
      alleles <- stringr::str_split(data_row, pattern = "', '")[[1]]  
    
      # Add vector to list
      temp_list[[row]] <- alleles
      }
    }
   
  # Set list item names to names of databases
  names(temp_list) <- databases
  
  # Remove NA entries
  temp_list2 <- temp_list[!is.na(temp_list)]
  
  # Add list for gene to overall list
  covered_alleles_list[[ gene ]] <- temp_list2
  
}

# Produce Venn diagrams
for (gene in genes) {
  
  # Get alleles covered by GeT-RM and Star Allele Search
  data_gene <- covered_alleles_list[[gene]]
  
  # Get all defined alleles from PharmVar
  all_alleles <- defined_alleles_list[[gene]]
  data_gene[[ "PharmVar" ]] <- all_alleles
  
  # Construct Venn diagram
  
  # If gene is not missing from any dataset
  if (length(data_gene) == 3) {
    venn.diagram(x = data_gene, main = gene, main.cex = 2.5,
                 category.names = c("Get-RM", "Star Allele Search", "PharmVar"),
                 filename = paste0('Plots/VennDiagram_', gene, '.png'),
                 cat.default.pos = "text", margin = 0.05, cat.cex = 1.6, cex = 2,
                 col=c("#440154ff", '#21908dff', '#fde725ff'), fontface = "bold", 
                 fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), 
                          alpha('#fde725ff',0.3)), imagetype = "png",
                 print.mode = c("raw"),
                 cat.just = list(c(0.3,0), c(0.6, 0), c(0.4,0)))
  } 
  # If gene is missing from GeT-RM
  else if (names(data_gene)[1] == "StarAlleleSearch") {
    venn.diagram(x = data_gene, main = gene, main.cex = 2.5,
                 category.names = c("Star Allele Search", "PharmVar"),
                 filename = paste0('Plots/VennDiagram_', gene, '.png'),
                 cat.default.pos = "text", margin = 0.05, cat.cex = 1.6, cex = 2,
                 col=c('#21908dff', '#fde725ff'), fontface = "bold", 
                 fill = c(alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
                 print.mode = c("raw"), imagetype = "png",
                 cat.just = list(c(0.6, 0), c(0.4,0)))
  }
  # If gene is missing from Star Allele Search
  else {
    venn.diagram(x = data_gene, main = gene, main.cex = 2.5,
                 category.names = c("Get-RM", "PharmVar"),
                 filename = paste0('Plots/VennDiagram_', gene, '.png'),
                 cat.default.pos = "text", margin = 0.05, cat.cex = 1.6, cex = 2,
                 col=c("#440154ff", '#fde725ff'), fontface = "bold", 
                 fill = c(alpha("#440154ff",0.3), alpha('#fde725ff',0.3)),
                 print.mode = c("raw"), imagetype = "png")
  }
}

################################
### DATABASES: MEAN COVERAGE ###
################################

### First: with equal weight given to all genes

apply(cov_databases, MARGIN = 2, FUN = mean)

# Also compute standard deviation!

apply(cov_databases, MARGIN = 2, FUN = sd)

### Next: using weights given by the total number of star alleles for each gene 

databases <- colnames(cov_databases)

# Find total number of defined star alleles for all genes
genes <- row.names(cov_databases)
total_num_alleles <- c(1:length(genes))
index <- 1
for (gene in genes) {
  alleles_gene <- defined_alleles_list[[ gene ]]
  total_num_alleles[index] <- length(alleles_gene)
  index <- index + 1
}

TotalAlleles <- data.frame("Gene" = genes, "NumAlleles" = total_num_alleles)

CoverageWeightedAverageDB <- function(database, include_missing_genes = FALSE) {
  # Compute weighted average of star allele coverage of all pharmacogenes in
  # PharmVar (excluding DPYD) for a specific database. An option is
  # provided to include genes that are not implemented by the caller in the
  # calculation as well.
  
  # Obtain allele coverages for all genes
  coverages <- cov_databases[, database]
  coverage_data <- data.frame("Gene" = row.names(cov_databases),
                              "Coverage" = coverages)
  
  # Add total number of alleles 
  coverage_data <- merge(coverage_data, TotalAlleles, by = "Gene")
  
  # If missing genes should not be included, remove them from the data frame
  if (!include_missing_genes) {
    coverage_data <- coverage_data %>% dplyr::filter(Coverage > 0)
  }
  
  # Compute numerator of weighted average
  numerator <- sum(coverage_data$Coverage * coverage_data$NumAlleles)
  
  # Compute denominator of weighted average
  denominator <- sum(coverage_data$NumAlleles)
  
  # Compute and return weighted average
  return(numerator/denominator)
}

CoverageWeightedAverageDB("GeT.RM", include_missing_genes = TRUE)
CoverageWeightedAverageDB("Star.allele.search", include_missing_genes = TRUE)

CoverageWeightedSdDB <- function(database, include_missing_genes = FALSE) {
  # Compute weighted standard deviation of star allele coverage of all pharmacogenes in
  # PharmVar (excluding DPYD) for a specific database. An option is
  # provided to include genes that are not present in the database in the
  # calculation as well.
  
  # Obtain allele coverages for all genes
  coverages <- cov_databases[, database]
  coverage_data <- data.frame("Gene" = row.names(cov_databases),
                              "Coverage" = coverages)
  
  # Add total number of alleles 
  coverage_data <- merge(coverage_data, TotalAlleles, by = "Gene")
  
  # If missing genes should not be included, remove them from the data frame
  if (!include_missing_genes) {
    coverage_data <- coverage_data %>% dplyr::filter(Coverage > 0)
  }
  
  # Compute weighted mean
  mean_w <- CoverageWeightedAverage(database, include_missing_genes)
  
  # Compute squared differences between weighted mean and coverages
  mean_w_diff <- (coverage_data$Coverage - mean_w)^2
  coverage_data$Squared_Diff <- mean_w_diff
  
  # Compute numerator of weighted sd
  numerator <- sum(coverage_data$Squared_Diff * coverage_data$NumAlleles)
  
  # Compute denominator of weighted sd
  denominator <- sum(coverage_data$NumAlleles) - 1
  
  # Compute and return weighted sd
  return(sqrt(numerator/denominator))
}

CoverageWeightedSdDB("GeT.RM", include_missing_genes = TRUE)
CoverageWeightedSdDB("Star.allele.search", include_missing_genes = TRUE)

##########################################
### STAR ALLELE CALLERS: MEAN COVERAGE ###
##########################################

### Compute weights given by the total number of star alleles for each gene 

callers <- c("PyPGx", "StellarPGx", "PharmCAT", "ursaPGx", "PAnno", "Aldy")

# Find total number of defined star alleles for all genes
genes <- row.names(cov_callers)
total_num_alleles <- c(1:length(genes))
index <- 1
for (gene in genes) {
  alleles_gene <- defined_alleles_list[[ gene ]]
  total_num_alleles[index] <- length(alleles_gene)
  index <- index + 1
}

TotalAlleles <- data.frame("Gene" = genes, "NumAlleles" = total_num_alleles)

CoverageWeightedAverageCaller <- function(caller, include_missing_genes = FALSE) {
  # Compute weighted average of star allele coverage of all pharmacogenes in
  # PharmVar (excluding DPYD) for a specific diplotype caller. An option is
  # provided to include genes that are not implemented by the caller in the
  # calculation as well.
  
  # Obtain allele coverages for all genes
  col_name <- paste0(caller, "_Coverage")
  coverages <- cov_callers[, col_name]
  coverage_data <- data.frame("Gene" = row.names(cov_callers),
                              "Coverage" = coverages)
  
  # Add total number of alleles 
  coverage_data <- merge(coverage_data, TotalAlleles, by = "Gene")
  
  # If missing genes should not be included, remove them from the data frame
  if (!include_missing_genes) {
    coverage_data <- coverage_data %>% dplyr::filter(Coverage > 0)
  }
  
  # Compute numerator of weighted average
  numerator <- sum(coverage_data$Coverage * coverage_data$NumAlleles)
  
  # Compute denominator of weighted average
  denominator <- sum(coverage_data$NumAlleles)
  
  # Compute and return weighted average
  return(numerator/denominator)
}

# Collect data required for plot
WA_callers_no_missing <- c()
WA_callers_missing <- c()

for (caller in callers) {
  cov1 <- CoverageWeightedAverageCaller(caller, include_missing_genes = FALSE)
  cov2 <- CoverageWeightedAverageCaller(caller, include_missing_genes = TRUE)
  WA_callers_no_missing <- c(WA_callers_no_missing, cov1)
  WA_callers_missing <- c(WA_callers_missing, cov2)
}

WA_callers_df <- data.frame(
  "Caller" = rep(callers, times = 2),
  "WeightedAverage" = c(WA_callers_missing, WA_callers_no_missing),
  "MissingInc" = factor(rep(c("Yes", "No"), each = length(callers)),
                        levels = c("Yes", "No"))
  )

labs_WA <- scales::percent(WA_callers_df$WeightedAverage, accuracy = .1)
WA_callers_df$Labels <- labs_WA

### Create plot

ggplot(data = WA_callers_df, aes(fill = MissingInc, color = MissingInc,
                                 y = WeightedAverage, x = Caller, label = Labels)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) + theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  geom_text(show.legend = FALSE, position = position_dodge(width = 0.8), 
            vjust = -0.5, size = 5.5, fontface = "bold") + 
  labs(y = "Weighted Average(*Allele Coverage)", fill = "Missing genes included", 
       color = "Missing genes included", x = "Star Allele Caller") +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold.italic"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 16, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size = 20)) +
  ggtitle("Weighted averages of allele coverage for six different diplotype callers") +
  scale_fill_manual(values = c("darkolivegreen3", "firebrick2")) +
  scale_color_manual(values = c("darkolivegreen4", "firebrick4"))

########################################
### STAR ALLELE CALLERS: UPSET PLOTS ###
########################################

# Read in data
callable_alleles <- read.csv("Callable_alleles.csv", sep = "\t")

# Get all genes to be included
genes <- c("CYP2B6", "CYP2A6", "CYP2D6", "CYP4F2", "NAT2")
callable_alleles_sub <- callable_alleles[, colnames(callable_alleles) %in% genes]

# Get caller names 
callers <- callable_alleles$Caller

# Initiate lists
callable_alleles_list <- list()

# Gather callable star alleles in a list of lists (one per gene)
for (gene in genes) {
  
  # Define list to store callable alleles for the gene as vectors
  temp_list <- list()
  
  for (row in 1:5) {
  
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

  # Add ursaPGx callable alleles
  temp_list[[6]] <- ursaPGx_callable_alleles_list[[ gene ]]
  
  # Set list item names to names of callers
  names(temp_list) <- c(callers, "ursaPGx")
  
  # Remove NA entries
  temp_list2 <- temp_list[!is.na(temp_list)]
  
  # Add list for gene to overall list
  callable_alleles_list[[ gene ]] <- temp_list2
  
}

## Make an upset plot for each of the genes
upset_plots <- list()
index <- 1
for (gene in genes) {
  
  # Obtain callable alleles for all callers
  data_gene <- callable_alleles_list[[gene]]
  
  # Obtain total number of defined star alleles for the gene (universal set)
  all_alleles <- defined_alleles_list[[gene]]

  # Preparation
  m2 <- list_to_matrix(data_gene, universal_set = all_alleles)
  m2 <- m2 == 1
  m2 <- as.data.frame(m2)
  gene_callers <- colnames(m2)
  
  # Create upset plot
  if (gene %in% c("CYP2B6", "CYP2D6")) { # No complement set for CYP2B6 and CYP2D6
    p2 <- upset(m2, gene_callers, width_ratio = 0.1, intersections = 'all',
                mode = "inclusive_union", max_degree = 2, name = paste("Caller combinations"),
                wrap = TRUE, stripes = c("lightgoldenrodyellow"), set_sizes = FALSE,
                themes=upset_default_themes(text=element_text(face = "bold",
                                                              size = 20)),
                base_annotations = list(
                  'Intersection size'=
                    intersection_size(text_mapping=aes(label=paste0(round(
                      !!get_size_mode('inclusive_union') / !!nrow(m2) * 100, 1), "%",
                      "\n", "(", !!get_size_mode("inclusive_union"), ")")),
                      mode = "inclusive_union", text=list(size = 6, fontface="bold"),
                      bar_number_threshold = 1, mapping = aes(fill = "bars_color")) +
                    ylab("Percentage coverage\n(absolute coverage)") +
                    scale_y_continuous(limits = c(0,nrow(m2) + 30),  
                                       labels = scales::percent_format(scale = 100/nrow(m2)),
                                       breaks=c(0, 20, 40, 60, 80, 100) / 100 * nrow(m2)) +
                    annotate(geom="text", x=Inf, y=Inf, label=paste("Total:", nrow(m2)),
                             vjust=1, hjust=1, size = 8, fontface = "bold") +
                    scale_fill_manual(values=c("bars_color"="darkolivegreen3"), guide="none")
                  ),
                matrix = (
                  intersection_matrix(
                    segment = geom_segment(linewidth = 1,
                                           colour = "darkolivegreen4"),
                    outline_color = list(
                      active = "darkolivegreen",
                      inactive = "grey"
                    )
                  )
                  + scale_color_manual(
                    values=c('TRUE'='darkolivegreen4', 'FALSE'='grey'),
                    guide = "none"
                    )
                  )
                ) + ggtitle(paste("Combined coverages of star allele callers for", gene)) +
      theme(plot.title = element_text(size = 25, face = "bold"))
  }
  else {
    p2 <- upset(m2, gene_callers, width_ratio = 0.1, intersections = 'all',
                mode = "inclusive_union", max_degree = 2, name = paste("Caller combinations"),
                wrap = TRUE, stripes = c("lightgoldenrodyellow"), set_sizes = FALSE,
                themes=upset_default_themes(text=element_text(face = "bold",
                                                              size = 20)),
                base_annotations = list(
                  'Intersection size'=
                    intersection_size(text_mapping=aes(label=paste0(round(
                      !!get_size_mode('inclusive_union') / !!nrow(m2) * 100, 1), "%",
                      "\n", "(", !!get_size_mode("inclusive_union"), ")")),
                      mode = "inclusive_union", text=list(size = 6, fontface="bold"),
                      bar_number_threshold = 1, mapping = aes(fill = "bars_color")) +
                    ylab("Percentage coverage\n(absolute coverage)") +
                    scale_y_continuous(limits = c(0,nrow(m2) + 1),  
                                       labels = scales::percent_format(scale = 100/nrow(m2)),
                                       breaks=c(0, 20, 40, 60, 80, 100) / 100 * nrow(m2)) +
                    annotate(geom="text", x=Inf, y=Inf, label=paste("Total:", nrow(m2)),
                             vjust=1, hjust=1, size = 8, fontface = "bold") +
                    scale_fill_manual(values=c("bars_color"="darkolivegreen3"), guide="none")
                  ),
                queries = list(
                  upset_query(
                    intersect = NA, color = "red3", fill = "red3"
                    )
                  ),
                matrix = (
                  intersection_matrix(
                    segment = geom_segment(linewidth = 1,
                                           colour = "darkolivegreen4"),
                    outline_color = list(
                      active = "darkolivegreen",
                      inactive = "grey"
                    )
                  )
                  + scale_color_manual(
                    values=c('TRUE'='darkolivegreen4', 'FALSE'='grey'),
                    guide = "none"
                    )
                  )
                ) + ggtitle(paste("Combined coverages of star allele callers for", gene)) +
      theme(plot.title = element_text(size = 25, face = "bold"))
  }
  upset_plots[[index]] <- p2
  index <- index + 1
}
