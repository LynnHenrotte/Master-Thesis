library("tidyverse")
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/Reading material")

data <- read.delim("relationships.tsv")

unique(data$Entity1_type)
unique(data$Entity2_type)

# Get gene-drug pairs from data
dataGeneDrug <- data %>% filter(Entity1_type == "Gene",
                                Entity2_type == "Chemical")

dataDrugGene <- data %>% filter(Entity1_type == "Chemical",
                                Entity2_type == "Gene")

# => same associations! (just entity 1 and 2 switched)

# Only select gene-drug pairs for which the evidence supports an association

dataAssociations <- dataGeneDrug %>% filter(Association == "associated")

# Get unique gene names
genes <- unique(dataAssociations$Entity1_name)
length(genes)

data2 <- read.delim("clinical_annotations.tsv")
unique(data2$Level.of.Evidence)

data2_sub <- data2 %>% filter(Level.of.Evidence %in% c("1A", "1B", '2A', '2B'))
length(unique(data2_sub$Gene))
