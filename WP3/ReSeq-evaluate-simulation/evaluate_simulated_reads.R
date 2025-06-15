#####################
### MASTER THESIS ###
#####################

### WP3: Simulation of Illumina short reads ###
### Lynn Henrotte ###

##----------------------------------------------------------------------------##

### LOAD PACKAGES ###

library("ggplot2")
library("tidyverse")
library("ggpubr")

##----------------------------------------------------------------------------##

##----------------------------------------------------------------------------##

###----- CYP2C9 -----###

##----------------------------------------------------------------------------##

### Per base average quality scores ###

cyp2c9_5x_qscores_files <- list.files("ReSeq_CYP2C9_5x_evaluation", 
                                       pattern = "*average_qscore_per_base.txt", full.names=TRUE)

cyp2c9_15x_qscores_files <- list.files("ReSeq_CYP2C9_15x_evaluation", 
                                       pattern = "*average_qscore_per_base.txt", full.names=TRUE)

cyp2c9_30x_qscores_files <- list.files("ReSeq_CYP2C9_30x_evaluation", 
                                       pattern = "*average_qscore_per_base.txt", full.names=TRUE)

### Manage and combine data ###

## For 5x coverage ##
for (allele in 1:length(cyp2c9_5x_qscores_files)) {
  
  qscore_file <- cyp2c9_5x_qscores_files[allele]
  qscore_data <- read.table(qscore_file, sep = "\t")
  
  read1_qscores <- qscore_data$V2 # change to V3 for log-transformed scores
  read2_qscores <- qscore_data$V5 # change to V6 for log-transformed scores
  positions <- qscore_data$V1
  
  name <- rep(paste0("CYP2C9*", allele), times = length(positions))
  label <- rep("Synthetic", times = length(positions))
  coverage <- rep("5x", times = length(positions))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_5x_qscores <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                                     "Read1_Qscore" = read1_qscores, "Read2_Qscore" = read2_qscores,
                                     "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                            "Read1_Qscore" = read1_qscores, "Read2_Qscore" = read2_qscores,
                            "Coverage" = coverage)

    cyp2c9_5x_qscores <- rbind(cyp2c9_5x_qscores, data_temp)
  }
}

# Add average quality score without making distinction between read 1 and 2
cyp2c9_5x_qscores <- cyp2c9_5x_qscores %>% mutate("Joint_Qscore" = (Read1_Qscore + Read2_Qscore)/2)

## For 15x coverage ##
for (allele in 1:length(cyp2c9_15x_qscores_files)) {
  
  qscore_file <- cyp2c9_15x_qscores_files[allele]
  qscore_data <- read.table(qscore_file, sep = "\t")
  
  read1_qscores <- qscore_data$V2 # change to V3 for log-transformed scores
  read2_qscores <- qscore_data$V5 # change to V6 for log-transformed scores
  positions <- qscore_data$V1
  
  name <- rep(paste0("CYP2C9*", allele), times = length(positions))
  label <- rep("Synthetic", times = length(positions))
  coverage <- rep("15x", times = length(positions))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_15x_qscores <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                                     "Read1_Qscore" = read1_qscores, "Read2_Qscore" = read2_qscores,
                                     "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                            "Read1_Qscore" = read1_qscores, "Read2_Qscore" = read2_qscores,
                            "Coverage" = coverage)
    
    cyp2c9_15x_qscores <- rbind(cyp2c9_15x_qscores, data_temp)
  }
}

# Add average quality score without making distinction between read 1 and 2
cyp2c9_15x_qscores <- cyp2c9_15x_qscores %>% mutate("Joint_Qscore" = (Read1_Qscore + Read2_Qscore)/2)

## For 30x coverage ##
for (allele in 1:length(cyp2c9_30x_qscores_files)) {
  
  qscore_file <- cyp2c9_30x_qscores_files[allele]
  qscore_data <- read.table(qscore_file, sep = "\t")
  
  read1_qscores <- qscore_data$V2 # change to V3 for log-transformed scores
  read2_qscores <- qscore_data$V5 # change to V6 for log-transformed scores
  positions <- qscore_data$V1
  
  name <- rep(paste0("CYP2C9*", allele), times = length(positions))
  label <- rep("Synthetic", times = length(positions))
  coverage <- rep("30x", times = length(positions))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_30x_qscores <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                                     "Read1_Qscore" = read1_qscores, "Read2_Qscore" = read2_qscores,
                                     "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                            "Read1_Qscore" = read1_qscores, "Read2_Qscore" = read2_qscores,
                            "Coverage" = coverage)
    
    cyp2c9_30x_qscores <- rbind(cyp2c9_30x_qscores, data_temp)
  }
}

# Add average quality score without making distinction between read 1 and 2
cyp2c9_30x_qscores <- cyp2c9_30x_qscores %>% mutate("Joint_Qscore" = (Read1_Qscore + Read2_Qscore)/2)

## Combine data from all coverages ##
cyp2c9_qscores <- rbind(cyp2c9_5x_qscores, cyp2c9_15x_qscores, cyp2c9_30x_qscores)

# Add unique identifier for each allele x coverage combination
cyp2c9_qscores <- cyp2c9_qscores %>% mutate("ID" = paste0(Allele, "_", Coverage))

# Make coverage a factor with ordered levels
cyp2c9_qscores$Coverage <- factor(cyp2c9_qscores$Coverage, levels = c("5x", "15x", "30x"))

## Plot average quality scores per position in the read (read 1 and 2 joined), for each allele

ggplot(data = cyp2c9_qscores, aes(x = Position, y = Joint_Qscore)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read", y = "Average quality score") +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

## Plot average quality scores per position in the first and second reads separately, for each allele

p_cyp2c9_qscores_read1 <- ggplot(data = cyp2c9_qscores, aes(x = Position, y = Read1_Qscore)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read 1", y = "Average quality score") + ylim(28,38) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

p_cyp2c9_qscores_read2 <- ggplot(data = cyp2c9_qscores, aes(x = Position, y = Read2_Qscore)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read 2", y = "Average quality score") + ylim(28, 38) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_cyp2c9_qscores_read1, p_cyp2c9_qscores_read2, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

## Take average over all simulated alleles (excluding CYP2C9*1)

cyp2c9_qscores_grouped <- cyp2c9_qscores %>% filter(Allele != "CYP2C9*1") %>%
  group_by(Position) %>% # Add Coverage to group_by to split up the lines by coverage
  summarize("Read1_Qscore" = mean(Read1_Qscore),
            "Read2_Qscore" = mean(Read2_Qscore),
            "Joint_Qscore" = mean(Joint_Qscore),
            .groups = "drop")

cyp2c9_qscores_grouped$Data <- "Synthetic"
cyp2c9_qscores_grouped$Group <- "Variants"

# Add CYP2C9*1 results

cyp2c9_1_qscores <- cyp2c9_qscores[cyp2c9_qscores$Allele == "CYP2C9*1",] %>%
  dplyr::select(!c(Allele, ID, Coverage)) %>% mutate(Group = "Wild type")

cyp2c9_1_qscores <- cyp2c9_1_qscores %>% group_by(Position) %>%
  summarize("Read1_Qscore" = mean(Read1_Qscore),
            "Read2_Qscore" = mean(Read2_Qscore),
            "Joint_Qscore" = mean(Joint_Qscore),
            .groups = "drop") %>%
  mutate("Data" = rep("Synthetic", times = 151),
         "Group" = rep("Wild type", times = 151))

cyp2c9_qscores_grouped <- rbind(cyp2c9_qscores_grouped, cyp2c9_1_qscores)

# Add real NovaSeq data results
cyp2c9_real_qscores_files <- list.files("evaluation_output", 
                                      pattern = "*average_qscore_per_base.txt", full.names=TRUE)

for (sample in 1:length(cyp2c9_real_qscores_files)) {
  
  qscore_file <- cyp2c9_real_qscores_files[sample]
  qscore_data <- read.table(qscore_file, sep = "\t")
  
  read1_qscores <- qscore_data$V2 # change to V3 for log-transformed scores
  read2_qscores <- qscore_data$V5 # change to V6 for log-transformed scores
  positions <- qscore_data$V1
  
  label <- rep("Real", times = length(positions))
  
  # Initiate data frame 
  if (sample == 1) {
    cyp2c9_real_qscores <- data.frame("Position" = positions, "Data" = label, "Read1_Qscore" = read1_qscores, 
                                      "Read2_Qscore" = read2_qscores)
  }
  else {
    data_temp <-  data.frame("Position" = positions, "Data" = label, "Read1_Qscore" = read1_qscores, 
                             "Read2_Qscore" = read2_qscores)
    
    cyp2c9_real_qscores <- rbind(cyp2c9_real_qscores, data_temp)
  }
}

# Add average quality score without making distinction between read 1 and 2
cyp2c9_real_qscores <- cyp2c9_real_qscores %>% mutate("Joint_Qscore" = (Read1_Qscore + Read2_Qscore)/2)

# Take average at each position over all samples
cyp2c9_real_qscores <- cyp2c9_real_qscores %>% group_by(Position) %>%
  summarize("Read1_Qscore" = mean(Read1_Qscore),
            "Read2_Qscore" = mean(Read2_Qscore),
            "Joint_Qscore" = mean(Joint_Qscore)) %>%
  mutate("Data" = "Real", "Group" = "Real")

# Add to simulated data
cyp2c9_qscores_grouped <- rbind(cyp2c9_qscores_grouped, cyp2c9_real_qscores)

## Plot per base quality scores, with average taken over all alleles with variants

p_qscores <- ggplot(data = cyp2c9_qscores_grouped, aes(x = Position, y = Joint_Qscore)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read", y = "Average quality score") + 
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

## Plot average quality scores per position in the first and second reads separately

p_cyp2c9_qscores_grouped_read1 <- ggplot(data = cyp2c9_qscores_grouped, aes(x = Position, y = Read1_Qscore)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 1", y = "Average quality score") + ylim(29, 38) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))


p_cyp2c9_qscores_grouped_read2 <- ggplot(data = cyp2c9_qscores_grouped, aes(x = Position, y = Read2_Qscore)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 2", y = "Average quality score") + ylim(29, 38) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))


ggarrange(p_cyp2c9_qscores_grouped_read1, p_cyp2c9_qscores_grouped_read2, 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

##----------------------------------------------------------------------------##

### Per base error rates ###

##----------------------------------------------------------------------------##

cyp2c9_5x_error_rate_files <- list.files("ReSeq_CYP2C9_5x_evaluation", 
                                      pattern = "*error_rates_per_base.txt", full.names=TRUE)

cyp2c9_15x_error_rate_files <- list.files("ReSeq_CYP2C9_15x_evaluation", 
                                       pattern = "*error_rates_per_base.txt", full.names=TRUE)

cyp2c9_30x_error_rate_files <- list.files("ReSeq_CYP2C9_30x_evaluation", 
                                       pattern = "*error_rates_per_base.txt", full.names=TRUE)

### Manage and combine data ###

## For 5x coverage ##

for (allele in 1:length(cyp2c9_5x_error_rate_files)) {
  
  error_file <- cyp2c9_5x_error_rate_files[allele]
  error_data <- read.table(error_file, sep = "\t", comment.char = "", header = TRUE)
  
  # Extract per base substitution error rates    
  read1_sub <- error_data$Sub1
  read2_sub <- error_data$Sub2
  
  # Extract per base deletion error rates
  read1_del <- error_data$Del1
  read2_del <- error_data$Del2
  
  # Extract per base insertion error rates
  read1_ins <- error_data$Ins1
  read2_ins <- error_data$Ins2
  
  # Compute per base indel error rates
  read1_indel <- read1_del + read1_ins
  read2_indel <- read2_del + read2_ins
  
  # Compute per base error rates
  read1_err <- read1_sub + read1_indel
  read2_err <- read2_sub + read2_indel
  
  # Obtain other relevant information
  positions <- error_data$X.BaseNum
  name <- rep(paste0("CYP2C9*", allele), times = length(positions))
  label <- rep("Synthetic", times = length(positions))
  coverage <- rep("5x", times = length(positions))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_5x_error_rate <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                                       "Read1_Sub" = read1_sub, "Read1_Del" = read1_del, "Read1_Ins" = read1_ins,
                                       "Read1_Indel" = read1_indel, "Read2_Sub" = read2_sub, "Read2_Del" = read2_del, 
                                       "Read2_Ins" = read2_ins, "Read2_Indel" = read2_indel, "Read1_Err" = read1_err,
                                       "Read2_Err" = read2_err, "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                            "Read1_Sub" = read1_sub, "Read1_Del" = read1_del, "Read1_Ins" = read1_ins,
                            "Read1_Indel" = read1_indel, "Read2_Sub" = read2_sub, "Read2_Del" = read2_del, 
                            "Read2_Ins" = read2_ins, "Read2_Indel" = read2_indel, "Read1_Err" = read1_err,
                            "Read2_Err" = read2_err, "Coverage" = coverage)
    
    cyp2c9_5x_error_rate <- rbind(cyp2c9_5x_error_rate, data_temp)
  }
}

# Add together error rates from read 1 and read 2 (take average)
cyp2c9_5x_error_rate <- cyp2c9_5x_error_rate %>% 
  mutate("Joint_Sub" = (Read1_Sub + Read2_Sub)/2,
         "Joint_Del" = (Read1_Del + Read2_Del)/2,
         "Joint_Ins" = (Read1_Ins + Read2_Ins)/2,
         "Joint_Indel" = (Read1_Indel + Read2_Indel)/2,
         "Joint_Err" = (Read1_Err + Read2_Err)/2)

## For 15x coverage ##

for (allele in 1:length(cyp2c9_15x_error_rate_files)) {
  
  error_file <- cyp2c9_15x_error_rate_files[allele]
  error_data <- read.table(error_file, sep = "\t", comment.char = "", header = TRUE)
  
  # Extract per base substitution error rates    
  read1_sub <- error_data$Sub1
  read2_sub <- error_data$Sub2
  
  # Extract per base deletion error rates
  read1_del <- error_data$Del1
  read2_del <- error_data$Del2
  
  # Extract per base insertion error rates
  read1_ins <- error_data$Ins1
  read2_ins <- error_data$Ins2
  
  # Compute per base indel error rates
  read1_indel <- read1_del + read1_ins
  read2_indel <- read2_del + read2_ins
  
  # Compute per base error rates
  read1_err <- read1_sub + read1_indel
  read2_err <- read2_sub + read2_indel
  
  # Obtain other relevant information
  positions <- error_data$X.BaseNum
  name <- rep(paste0("CYP2C9*", allele), times = length(positions))
  label <- rep("Synthetic", times = length(positions))
  coverage <- rep("15x", times = length(positions))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_15x_error_rate <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                                        "Read1_Sub" = read1_sub, "Read1_Del" = read1_del, "Read1_Ins" = read1_ins,
                                        "Read1_Indel" = read1_indel, "Read2_Sub" = read2_sub, "Read2_Del" = read2_del, 
                                        "Read2_Ins" = read2_ins, "Read2_Indel" = read2_indel, "Read1_Err" = read1_err,
                                        "Read2_Err" = read2_err, "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                            "Read1_Sub" = read1_sub, "Read1_Del" = read1_del, "Read1_Ins" = read1_ins,
                            "Read1_Indel" = read1_indel, "Read2_Sub" = read2_sub, "Read2_Del" = read2_del, 
                            "Read2_Ins" = read2_ins, "Read2_Indel" = read2_indel, "Read1_Err" = read1_err,
                            "Read2_Err" = read2_err, "Coverage" = coverage)
    
    cyp2c9_15x_error_rate <- rbind(cyp2c9_15x_error_rate, data_temp)
  }
}

# Add together error rates from read 1 and read 2 (take average)
cyp2c9_15x_error_rate <- cyp2c9_15x_error_rate %>% 
  mutate("Joint_Sub" = (Read1_Sub + Read2_Sub)/2,
         "Joint_Del" = (Read1_Del + Read2_Del)/2,
         "Joint_Ins" = (Read1_Ins + Read2_Ins)/2,
         "Joint_Indel" = (Read1_Indel + Read2_Indel)/2,
         "Joint_Err" = (Read1_Err + Read2_Err)/2)

## For 30x coverage ##

for (allele in 1:length(cyp2c9_30x_error_rate_files)) {
  
  error_file <- cyp2c9_30x_error_rate_files[allele]
  error_data <- read.table(error_file, sep = "\t", comment.char = "", header = TRUE)
  
  # Extract per base substitution error rates    
  read1_sub <- error_data$Sub1
  read2_sub <- error_data$Sub2
  
  # Extract per base deletion error rates
  read1_del <- error_data$Del1
  read2_del <- error_data$Del2
  
  # Extract per base insertion error rates
  read1_ins <- error_data$Ins1
  read2_ins <- error_data$Ins2
  
  # Compute per base indel error rates
  read1_indel <- read1_del + read1_ins
  read2_indel <- read2_del + read2_ins
  
  # Compute per base error rates
  read1_err <- read1_sub + read1_indel
  read2_err <- read2_sub + read2_indel
  
  # Obtain other relevant information
  positions <- error_data$X.BaseNum
  name <- rep(paste0("CYP2C9*", allele), times = length(positions))
  label <- rep("Synthetic", times = length(positions))
  coverage <- rep("30x", times = length(positions))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_30x_error_rate <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                                       "Read1_Sub" = read1_sub, "Read1_Del" = read1_del, "Read1_Ins" = read1_ins,
                                       "Read1_Indel" = read1_indel, "Read2_Sub" = read2_sub, "Read2_Del" = read2_del, 
                                       "Read2_Ins" = read2_ins, "Read2_Indel" = read2_indel, "Read1_Err" = read1_err,
                                       "Read2_Err" = read2_err, "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("Position" = positions, "Data" = label, "Allele" = name, 
                            "Read1_Sub" = read1_sub, "Read1_Del" = read1_del, "Read1_Ins" = read1_ins,
                            "Read1_Indel" = read1_indel, "Read2_Sub" = read2_sub, "Read2_Del" = read2_del, 
                            "Read2_Ins" = read2_ins, "Read2_Indel" = read2_indel, "Read1_Err" = read1_err,
                            "Read2_Err" = read2_err, "Coverage" = coverage)
    
    cyp2c9_30x_error_rate <- rbind(cyp2c9_30x_error_rate, data_temp)
  }
}

# Add together error rates from read 1 and read 2 (take average)
cyp2c9_30x_error_rate <- cyp2c9_30x_error_rate %>% 
  mutate("Joint_Sub" = (Read1_Sub + Read2_Sub)/2,
         "Joint_Del" = (Read1_Del + Read2_Del)/2,
         "Joint_Ins" = (Read1_Ins + Read2_Ins)/2,
         "Joint_Indel" = (Read1_Indel + Read2_Indel)/2,
         "Joint_Err" = (Read1_Err + Read2_Err)/2)

## Combine data from all coverages ##
cyp2c9_error_rate <- rbind(cyp2c9_5x_error_rate, cyp2c9_15x_error_rate, cyp2c9_30x_error_rate)

# Add unique identifier for each allele x coverage combination
cyp2c9_error_rate <- cyp2c9_error_rate %>% mutate("ID" = paste0(Allele, "_", Coverage))

# Make coverage a factor with ordered levels
cyp2c9_error_rate$Coverage <- factor(cyp2c9_error_rate$Coverage, levels = c("5x", "15x", "30x"))

## Plot error rate per position in the read (read 1 and 2 joined), for each allele

# Overall error rate

ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Joint_Err)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read", y = "Overall error rate") + ylim(0, 0.06) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))


# Substitution error rate

ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Joint_Sub)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read", y = "Substitution error rate") + ylim(0, 0.06) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))


# Indel error rate

ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Joint_Indel)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read", y = "Indel error rate") + ylim(0,0.0025) +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold.italic"),
        legend.text = element_text(size = 16, face = "italic"),
        legend.title = element_text(size = 16, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size = 25))

## Plot average error rates per position in the first and second reads separately, for each allele

# Overall error rate

p_cyp2c9_err_read1 <- ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Read1_Err)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read 1", y = "Overall error rate") + ylim(0, 0.07) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))


p_cyp2c9_err_read2 <- ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Read2_Err)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read 2", y = "Overall error rate") + ylim(0, 0.07) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_cyp2c9_err_read1, p_cyp2c9_err_read2, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

# Substitution error rate

p_cyp2c9_err_sub_read1 <- ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Read1_Sub)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read 1", y = "Substitution error rate") + ylim(0, 0.07) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

p_cyp2c9_err_sub_read2 <- ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Read2_Sub)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read 2", y = "Substitution error rate") + ylim(0, 0.07) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_cyp2c9_err_sub_read1, p_cyp2c9_err_sub_read2, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

# Indel error rate

p_cyp2c9_err_indel_read1 <- ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Read1_Indel)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read 1", y = "Indel error rate") + ylim(0, 0.0025) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

p_cyp2c9_err_indel_read2 <- ggplot(data = cyp2c9_error_rate, aes(x = Position, y = Read2_Indel)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read 2", y = "Indel error rate") + ylim(0, 0.0025) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_cyp2c9_err_indel_read1, p_cyp2c9_err_indel_read2, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

## Take average over all simulated alleles (excluding CYP2C9*1)

cyp2c9_error_rate_grouped <- cyp2c9_error_rate %>% filter(Allele != "CYP2C9*1") %>%
  group_by(Position) %>% 
  summarize("Read1_Err" = mean(Read1_Err), "Read2_Err" = mean(Read2_Err),
            "Read1_Sub" = mean(Read1_Sub), "Read2_Sub" = mean(Read2_Sub),
            "Read1_Del" = mean(Read1_Del), "Read2_Del" = mean(Read2_Del),
            "Read1_Ins" = mean(Read1_Ins), "Read2_Ins" = mean(Read2_Ins),
            "Read1_Indel" = mean(Read1_Indel), "Read2_Indel" = mean(Read2_Indel),
            "Joint_Err" = mean(Joint_Err), "Joint_Sub" = mean(Joint_Sub),
            "Joint_Del" = mean(Joint_Del), "Joint_Ins" = mean(Joint_Ins),
            "Joint_Indel" = mean(Joint_Indel), .groups = "drop")

cyp2c9_error_rate_grouped$Data <- "Synthetic"
cyp2c9_error_rate_grouped$Group <- "Variants"

# Add CYP2C9*1 results

cyp2c9_1_error_rate <- cyp2c9_error_rate[cyp2c9_error_rate$Allele == "CYP2C9*1",] %>%
  dplyr::select(!c(Allele, Coverage, ID)) %>% mutate(Group = "Wild type")

cyp2c9_1_error_rate <- cyp2c9_1_error_rate %>% group_by(Position) %>%
  summarize("Read1_Err" = mean(Read1_Err), "Read2_Err" = mean(Read2_Err),
            "Read1_Sub" = mean(Read1_Sub), "Read2_Sub" = mean(Read2_Sub),
            "Read1_Del" = mean(Read1_Del), "Read2_Del" = mean(Read2_Del),
            "Read1_Ins" = mean(Read1_Ins), "Read2_Ins" = mean(Read2_Ins),
            "Read1_Indel" = mean(Read1_Indel), "Read2_Indel" = mean(Read2_Indel),
            "Joint_Err" = mean(Joint_Err), "Joint_Sub" = mean(Joint_Sub),
            "Joint_Del" = mean(Joint_Del), "Joint_Ins" = mean(Joint_Ins),
            "Joint_Indel" = mean(Joint_Indel), .groups = "drop") %>%
  mutate("Data" = rep("Synthetic", times = 151),
         "Group" = rep("Wild type", times = 151))

cyp2c9_error_rate_grouped <- rbind(cyp2c9_error_rate_grouped, cyp2c9_1_error_rate)

# Add real NovaSeq data results
cyp2c9_real_error_files <- list.files("evaluation_output", 
                                        pattern = "*error_rates_per_base.txt", full.names=TRUE)

for (sample in 1:length(cyp2c9_real_error_files)) {
  
  error_file <- cyp2c9_real_error_files[sample]
  error_data <- read.table(error_file, sep = "\t", comment.char = "", header = TRUE)
  
  # Extract per base substitution error rates    
  read1_sub <- error_data$Sub1
  read2_sub <- error_data$Sub2
  
  # Extract per base deletion error rates
  read1_del <- error_data$Del1
  read2_del <- error_data$Del2
  
  # Extract per base insertion error rates
  read1_ins <- error_data$Ins1
  read2_ins <- error_data$Ins2
  
  # Compute per base indel error rates
  read1_indel <- read1_del + read1_ins
  read2_indel <- read2_del + read2_ins
  
  # Compute per base error rates
  read1_err <- read1_sub + read1_indel
  read2_err <- read2_sub + read2_indel
  
  # Obtain other relevant information
  positions <- error_data$X.BaseNum
  label <- rep("Real", times = length(positions))

  # Initiate data frame 
  if (sample == 1) {
    cyp2c9_real_errors <- data.frame("Position" = positions, "Data" = label,
                                      "Read1_Sub" = read1_sub, "Read1_Del" = read1_del, "Read1_Ins" = read1_ins,
                                      "Read1_Indel" = read1_indel, "Read2_Sub" = read2_sub, "Read2_Del" = read2_del, 
                                      "Read2_Ins" = read2_ins, "Read2_Indel" = read2_indel, "Read1_Err" = read1_err,
                                      "Read2_Err" = read2_err)
  }
  else {
    data_temp <- data.frame("Position" = positions, "Data" = label,
                              "Read1_Sub" = read1_sub, "Read1_Del" = read1_del, "Read1_Ins" = read1_ins,
                              "Read1_Indel" = read1_indel, "Read2_Sub" = read2_sub, "Read2_Del" = read2_del, 
                              "Read2_Ins" = read2_ins, "Read2_Indel" = read2_indel, "Read1_Err" = read1_err,
                              "Read2_Err" = read2_err)
    
    cyp2c9_real_errors <- rbind(cyp2c9_real_errors, data_temp)
  }
}

# Add error rates without making distinction between read 1 and 2
cyp2c9_real_errors <- cyp2c9_real_errors %>% 
  mutate("Joint_Sub" = (Read1_Sub + Read2_Sub)/2,
         "Joint_Del" = (Read1_Del + Read2_Del)/2,
         "Joint_Ins" = (Read1_Ins + Read2_Ins)/2,
         "Joint_Indel" = (Read1_Indel + Read2_Indel)/2,
         "Joint_Err" = (Read1_Err + Read2_Err)/2)

# Take average at each position over all samples
cyp2c9_real_errors <- cyp2c9_real_errors %>% group_by(Position) %>%
  summarize("Read1_Err" = mean(Read1_Err), "Read2_Err" = mean(Read2_Err),
            "Read1_Sub" = mean(Read1_Sub), "Read2_Sub" = mean(Read2_Sub),
            "Read1_Del" = mean(Read1_Del), "Read2_Del" = mean(Read2_Del),
            "Read1_Ins" = mean(Read1_Ins), "Read2_Ins" = mean(Read2_Ins),
            "Read1_Indel" = mean(Read1_Indel), "Read2_Indel" = mean(Read2_Indel),
            "Joint_Err" = mean(Joint_Err), "Joint_Sub" = mean(Joint_Sub),
            "Joint_Del" = mean(Joint_Del), "Joint_Ins" = mean(Joint_Ins),
            "Joint_Indel" = mean(Joint_Indel), .groups = "drop") %>%
  mutate("Data" = "Real", "Group" = "Real")

# Add to simulated data
cyp2c9_error_rate_grouped <- rbind(cyp2c9_error_rate_grouped, cyp2c9_real_errors)

## Plot per base overall error rate, with average taken over all alleles with variants

# Overall error rate

p_error <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Joint_Err)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read", y = "Overall error rate") + ylim(0, 0.06) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))


# Substitution error rate

p_sub <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Joint_Sub)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read", y = "Substitution error rate") + ylim(0, 0.06) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

# Indel error rate

p_indel <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Joint_Indel)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read", y = "Indel error rate") + ylim(0, 0.002) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

# Indel error rate, ignoring wild type

cyp2c9_error_rate_grouped_no_wild_type <- cyp2c9_error_rate_grouped %>%
  filter(Group != "Wild type")

p_indel2 <- ggplot(data = cyp2c9_error_rate_grouped_no_wild_type, aes(x = Position, y = Joint_Indel)) + 
  geom_line(aes(group = Data, color = Data), linewidth = 1.5, alpha = 0.7) + theme_bw() +
  labs(x = "Position in read", y = "Indel error rate") + ylim(0, 0.002) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_sub, p_indel, nrow = 1, ncol = 2, common.legend = TRUE,
          legend = "bottom")

## Plot error rate per position in the first and second reads separately

# Overall error rate

p_cyp2c9_err_grouped_read1 <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Read1_Err)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 1", y = "Overall error rate") + ylim(0,0.065) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

p_cyp2c9_err_grouped_read2 <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Read2_Err)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 2", y = "Overall error rate") + ylim(0,0.065) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_cyp2c9_err_grouped_read1, p_cyp2c9_err_grouped_read2, 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

# Substitution error rate

p_cyp2c9_err_sub_grouped_read1 <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Read1_Sub)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 1", y = "Substitution error rate") + ylim(0,0.065) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

p_cyp2c9_err_sub_grouped_read2 <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Read2_Sub)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 2", y = "Substitution error rate") + ylim(0,0.065) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_cyp2c9_err_sub_grouped_read1, p_cyp2c9_err_sub_grouped_read2, 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

# Indel error rate

p_cyp2c9_err_indel_grouped_read1 <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Read1_Indel)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 1", y = "Indel error rate") + ylim(0,0.002) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

p_cyp2c9_err_indel_grouped_read2 <- ggplot(data = cyp2c9_error_rate_grouped, aes(x = Position, y = Read2_Indel)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 2", y = "Indel error rate") + ylim(0,0.002) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_cyp2c9_err_indel_grouped_read1, p_cyp2c9_err_indel_grouped_read2, 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

# Indel error rate, ignoring wild type 

cyp2c9_error_rate_grouped_no_wild_type <- cyp2c9_error_rate_grouped %>%
  filter(Group != "Wild type")

p_cyp2c9_err_indel_grouped_no_wildtype_read1 <- ggplot(data = cyp2c9_error_rate_grouped_no_wild_type, 
                                           aes(x = Position, y = Read1_Indel)) + 
  geom_line(aes(group = Data, color = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 1", y = "Indel error rate") + ylim(0,0.002) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

p_cyp2c9_err_indel_grouped_no_wildtype_read2 <- ggplot(data = cyp2c9_error_rate_grouped_no_wild_type, 
                                           aes(x = Position, y = Read2_Indel)) + 
  geom_line(aes(group = Data, color = Data), linewidth = 1.5, alpha = 0.8) + theme_bw() +
  labs(x = "Position in read 2", y = "Indel error rate") + ylim(0,0.002) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

ggarrange(p_cyp2c9_err_indel_grouped_no_wildtype_read1, 
          p_cyp2c9_err_indel_grouped_no_wildtype_read2, 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

##----------------------------------------------------------------------------##

### GC content ###

##----------------------------------------------------------------------------##

cyp2c9_5x_GC_files <- list.files("ReSeq_CYP2C9_5x_evaluation", 
                                 pattern = "*GC_content.txt", full.names=TRUE)

cyp2c9_15x_GC_files <- list.files("ReSeq_CYP2C9_15x_evaluation", 
                                  pattern = "*GC_content.txt", full.names=TRUE)

cyp2c9_30x_GC_files <- list.files("ReSeq_CYP2C9_30x_evaluation",
                                  pattern = "*GC_content.txt", full.names=TRUE)

### Manage and combine data ###

## For 5x coverage ##

for (allele in 1:length(cyp2c9_5x_GC_files)) {
  
  GC_file <- cyp2c9_5x_GC_files[allele]
  GC_data <- read.table(GC_file, sep = "\t", col.names = c("GC", "Count"))
  
  # Obtain number of read pairs for each percentage of GC content
  GC_count <- GC_data$Count
  GC_content <- GC_data$GC
  
  # Convert counts to percentages
  GC_percent <- GC_count / sum(GC_count)
  
  # Obtain other relevant information
  name <- rep(paste0("CYP2C9*", allele), times = length(GC_percent))
  label <- rep("Synthetic", times = length(GC_percent))
  coverage <- rep("5x", times = length(GC_percent))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_5x_GC <- data.frame("GC_content" = GC_content, "Count" = GC_count, 
                               "Percent" = GC_percent, "Data" = label, "Allele" = name, 
                               "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("GC_content" = GC_content, "Count" = GC_count, 
                            "Percent" = GC_percent, "Data" = label, "Allele" = name, 
                            "Coverage" = coverage)
    
    cyp2c9_5x_GC <- rbind(cyp2c9_5x_GC, data_temp)
  }
}

## For 15x coverage ##

for (allele in 1:length(cyp2c9_15x_GC_files)) {
  
  GC_file <- cyp2c9_15x_GC_files[allele]
  GC_data <- read.table(GC_file, sep = "\t", col.names = c("GC", "Count"))
  
  # Obtain number of read pairs for each percentage of GC content
  GC_count <- GC_data$Count
  GC_content <- GC_data$GC
  
  # Convert counts to percentages
  GC_percent <- GC_count / sum(GC_count)
  
  # Obtain other relevant information
  name <- rep(paste0("CYP2C9*", allele), times = length(GC_percent))
  label <- rep("Synthetic", times = length(GC_percent))
  coverage <- rep("15x", times = length(GC_percent))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_15x_GC <- data.frame("GC_content" = GC_content, "Count" = GC_count, 
                               "Percent" = GC_percent, "Data" = label, "Allele" = name, 
                               "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("GC_content" = GC_content, "Count" = GC_count, 
                            "Percent" = GC_percent, "Data" = label, "Allele" = name, 
                            "Coverage" = coverage)
    
    cyp2c9_15x_GC <- rbind(cyp2c9_15x_GC, data_temp)
  }
}

## For 30x coverage ##

for (allele in 1:length(cyp2c9_30x_GC_files)) {
  
  GC_file <- cyp2c9_30x_GC_files[allele]
  GC_data <- read.table(GC_file, sep = "\t", col.names = c("GC", "Count"))
  
  # Obtain number of read pairs for each percentage of GC content
  GC_count <- GC_data$Count
  GC_content <- GC_data$GC
  
  # Convert counts to percentages
  GC_percent <- GC_count / sum(GC_count)
  
  # Obtain other relevant information
  name <- rep(paste0("CYP2C9*", allele), times = length(GC_percent))
  label <- rep("Synthetic", times = length(GC_percent))
  coverage <- rep("30x", times = length(GC_percent))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_30x_GC <- data.frame("GC_content" = GC_content, "Count" = GC_count, 
                               "Percent" = GC_percent, "Data" = label, "Allele" = name, 
                               "Coverage" = coverage)
  }
  else {
    data_temp <- data.frame("GC_content" = GC_content, "Count" = GC_count, 
                            "Percent" = GC_percent, "Data" = label, "Allele" = name, 
                            "Coverage" = coverage)
    
    cyp2c9_30x_GC <- rbind(cyp2c9_30x_GC, data_temp)
  }
}

## Combine data from all coverages ##
cyp2c9_GC <- rbind(cyp2c9_5x_GC, cyp2c9_15x_GC, cyp2c9_30x_GC)

# Add unique identifier for each allele x coverage combination
cyp2c9_GC <- cyp2c9_GC %>% mutate("ID" = paste0(Allele, "_", Coverage))

# Make coverage a factor with ordered levels
cyp2c9_GC$Coverage <- factor(cyp2c9_GC$Coverage, levels = c("5x", "15x", "30x"))

## Plot per sequence GC content distribution for each allele ##

ggplot(data = cyp2c9_GC, aes(x = GC_content, y = Percent)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + 
  theme_bw() +
  labs(x = "Mean GC content (%)", y = "Frequency") + ylim(0, 0.085) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

## Take average over all simulated alleles (excluding CYP2C9*1)

cyp2c9_GC_grouped <- cyp2c9_GC %>% 
  filter(Allele != "CYP2C9*1") %>% group_by(GC_content) %>% 
  summarize("Percent" = mean(Percent), .groups = "drop")

cyp2c9_GC_grouped$Data <- "Synthetic"
cyp2c9_GC_grouped$Group <- "Variants"

# Add CYP2C9*1 results

cyp2c9_1_GC <- cyp2c9_GC[cyp2c9_GC$Allele == "CYP2C9*1",] %>%
  group_by(GC_content) %>%
  summarize("Percent" = mean(Percent), .groups = "drop") %>%
  mutate(Group = "Wild type", Data = "Synthetic")
cyp2c9_GC_grouped <- rbind(cyp2c9_GC_grouped, cyp2c9_1_GC)

# Add real NovaSeq data results
cyp2c9_real_GC_files <- list.files("evaluation_output", pattern = "*GC_content.txt", full.names=TRUE)

for (sample in 1:length(cyp2c9_real_GC_files)) {
  
  GC_file <- cyp2c9_real_GC_files[sample]
  GC_data <- read.table(GC_file, sep = "\t", col.names = c("GC", "Count"))
  
  # Obtain number of read pairs for each percentage of GC content
  GC_count <- GC_data$Count
  GC_content <- GC_data$GC
  
  # Convert counts to percentages
  GC_percent <- GC_count / sum(GC_count)
  
  # Define label
  label <- rep("Real", times = length(GC_percent))
  
  # Initiate data frame 
  if (sample == 1) {
    cyp2c9_real_GC_all <- data.frame("GC_content" = GC_content, "Count" = GC_count, 
                                "Percent" = GC_percent, "Data" = label)
  }
  else {
    data_temp <- data.frame("GC_content" = GC_content, "Count" = GC_count, 
                            "Percent" = GC_percent, "Data" = label)
    
    cyp2c9_real_GC_all <- rbind(cyp2c9_real_GC_all, data_temp)
  }
}

# Take average of frequencies over all samples 

cyp2c9_real_GC <- cyp2c9_real_GC_all %>% group_by(GC_content) %>% 
  summarize("Percent" = mean(Percent), .groups = "drop") %>%
  mutate(Data = "Real", Group = "Real")

# Add real data results to simulated data results

cyp2c9_GC_grouped <- rbind(cyp2c9_GC_grouped, cyp2c9_real_GC)

## Plot per sequence GC content distribution for wild type and variants (average) ##

ggplot(data = cyp2c9_GC_grouped, aes(x = GC_content, y = Percent)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), 
            linewidth = 1.5, alpha = 0.7) + 
  theme_bw() + labs(x = "Mean GC content (%)", y = "Frequency") + ylim(0, 0.085) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

## Take average over all simulated alleles, including CYP2C9*1

cyp2c9_GC_grouped2 <- cyp2c9_GC %>% group_by(GC_content) %>% 
  summarize("Percent" = mean(Percent), .groups = "drop")

cyp2c9_GC_grouped2$Data <- "Synthetic"

## Add real NovaSeq data results

cyp2c9_GC_grouped2 <- rbind(cyp2c9_GC_grouped2, cyp2c9_real_GC %>% dplyr::select(!c(Group)))

## Plot per sequence GC content distribution for simulated data and real data ##

p_GC <- ggplot(data = cyp2c9_GC_grouped2, aes(x = GC_content, y = Percent)) + 
  geom_line(aes(group = Data, color = Data), 
            linewidth = 1.5, alpha = 0.7) + 
  theme_bw() + labs(x = "Mean GC content (%)", y = "Frequency") + ylim(0, 0.085) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

# Make vectors expanding GC content for real and simulated data
cyp2c9_GC_total_counts <- rbind(cyp2c9_GC %>% dplyr::select(c(GC_content, Count, Data)),
                                cyp2c9_real_GC_all %>% dplyr::select(c(GC_content, Count, Data)))
                                
cyp2c9_GC_total_counts <- cyp2c9_GC_total_counts %>% group_by(Data, GC_content) %>%
  summarize("Total-Count" = sum(Count), .groups = "drop")

GC_content_synthetic_reads <- c()
GC_content_real_reads <- c()

for (i in 0:100) {
  
  # Get real data count
  count_real <- (cyp2c9_GC_total_counts %>% 
                   filter(Data == "Real", GC_content == i))$`Total-Count`
  
  # Get simulated data count
  count_synthetic <- (cyp2c9_GC_total_counts %>% 
                        filter(Data == "Synthetic", GC_content == i))$`Total-Count`
  
  # Add to vectors
  GC_content_real_reads <- c(GC_content_real_reads, rep(i, times = count_real))
  GC_content_synthetic_reads <- c(GC_content_synthetic_reads, rep(i, times = count_synthetic))
  
}

# Compute means directly
mean(GC_content_real_reads)
mean(GC_content_synthetic_reads)

# Compute variances directly
GC_var_real <- var(GC_content_real_reads)
GC_var_synthetic <- var(GC_content_synthetic_reads)

sqrt(GC_var_real)
sqrt(GC_var_synthetic)

##----------------------------------------------------------------------------##

### Insert size distribution ###

##----------------------------------------------------------------------------##

cyp2c9_5x_insert_files <- list.files("ReSeq_CYP2C9_5x_evaluation", 
                                     pattern = "*insert_sizes.txt", full.names=TRUE)

cyp2c9_15x_insert_files <- list.files("ReSeq_CYP2C9_15x_evaluation", 
                                      pattern = "*insert_sizes.txt", full.names=TRUE)

cyp2c9_30x_insert_files <- list.files("ReSeq_CYP2C9_30x_evaluation",
                                      pattern = "*insert_sizes.txt", full.names=TRUE)

cyp2c9_real_insert_files <- list.files("evaluation_output", 
                                       pattern = "*insert_sizes.txt", full.names=TRUE)

### Find all insert sizes

insert_sizes_all <- c()

for (file in c(cyp2c9_5x_insert_files, cyp2c9_15x_insert_files, cyp2c9_30x_insert_files,
               cyp2c9_real_insert_files)) {
  
  data <- read.table(file, sep = "\t", col.names = c("Insert_Size", "Count"))
  insert_sizes_data <- data$Insert_Size
  
  insert_sizes_all <- union(insert_sizes_all, insert_sizes_data)
}

IS_data_join <- data.frame("Insert_Size" = insert_sizes_all) %>% arrange(Insert_Size)

### Manage and combine data ###

## For 5x coverage ##

for (allele in 1:length(cyp2c9_5x_insert_files)) {
  
  insert_file <- cyp2c9_5x_insert_files[allele]
  insert_data <- read.table(insert_file, sep = "\t", col.names = c("Insert_Size", "Count"))
  
  # Obtain counts and percentages of each insert size
  IS_count <- insert_data$Count
  insert_data$Percent <- IS_count / sum(IS_count)
  insert_sizes <- insert_data$Insert_Size
  
  # Obtain other relevant information
  name <- rep(paste0("CYP2C9*", allele), times = length(insert_sizes_all))
  label <- rep("Synthetic", times = length(insert_sizes_all))
  coverage <- rep("5x", times = length(insert_sizes_all))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_5x_insert <- merge(IS_data_join, insert_data, by = "Insert_Size", all.x = TRUE)
    cyp2c9_5x_insert$Coverage <- coverage
    cyp2c9_5x_insert$Data <- label
    cyp2c9_5x_insert$Allele <- name
    
  }
  else {
    data_temp <- merge(IS_data_join, insert_data, by = "Insert_Size", all.x = TRUE)
    data_temp$Coverage <- coverage
    data_temp$Data <- label
    data_temp$Allele <- name
    
    cyp2c9_5x_insert <- rbind(cyp2c9_5x_insert, data_temp)
  }
}

# Set all NA values to zero
cyp2c9_5x_insert[is.na(cyp2c9_5x_insert)] <- 0

## For 15x coverage ##

for (allele in 1:length(cyp2c9_15x_insert_files)) {
  
  insert_file <- cyp2c9_15x_insert_files[allele]
  insert_data <- read.table(insert_file, sep = "\t", col.names = c("Insert_Size", "Count"))
  
  # Obtain counts and percentages of each insert size
  IS_count <- insert_data$Count
  insert_data$Percent <- IS_count / sum(IS_count)
  insert_sizes <- insert_data$Insert_Size
  
  # Obtain other relevant information
  name <- rep(paste0("CYP2C9*", allele), times = length(insert_sizes_all))
  label <- rep("Synthetic", times = length(insert_sizes_all))
  coverage <- rep("15x", times = length(insert_sizes_all))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_15x_insert <- merge(IS_data_join, insert_data, by = "Insert_Size", all.x = TRUE)
    cyp2c9_15x_insert$Coverage <- coverage
    cyp2c9_15x_insert$Data <- label
    cyp2c9_15x_insert$Allele <- name
    
  }
  else {
    data_temp <- merge(IS_data_join, insert_data, by = "Insert_Size", all.x = TRUE)
    data_temp$Coverage <- coverage
    data_temp$Data <- label
    data_temp$Allele <- name
    
    cyp2c9_15x_insert <- rbind(cyp2c9_15x_insert, data_temp)
  }
}

# Set all NA values to zero
cyp2c9_15x_insert[is.na(cyp2c9_15x_insert)] <- 0

## For 30x coverage ##

for (allele in 1:length(cyp2c9_30x_insert_files)) {
  
  insert_file <- cyp2c9_30x_insert_files[allele]
  insert_data <- read.table(insert_file, sep = "\t", col.names = c("Insert_Size", "Count"))
  
  # Obtain counts and percentages of each insert size
  IS_count <- insert_data$Count
  insert_data$Percent <- IS_count / sum(IS_count)
  insert_sizes <- insert_data$Insert_Size
  
  # Obtain other relevant information
  name <- rep(paste0("CYP2C9*", allele), times = length(insert_sizes_all))
  label <- rep("Synthetic", times = length(insert_sizes_all))
  coverage <- rep("30x", times = length(insert_sizes_all))
  
  # Initiate data frame 
  if (allele == 1) {
    cyp2c9_30x_insert <- merge(IS_data_join, insert_data, by = "Insert_Size", all.x = TRUE)
    cyp2c9_30x_insert$Coverage <- coverage
    cyp2c9_30x_insert$Data <- label
    cyp2c9_30x_insert$Allele <- name
    
  }
  else {
    data_temp <- merge(IS_data_join, insert_data, by = "Insert_Size", all.x = TRUE)
    data_temp$Coverage <- coverage
    data_temp$Data <- label
    data_temp$Allele <- name
    
    cyp2c9_30x_insert <- rbind(cyp2c9_30x_insert, data_temp)
  }
}

# Set all NA values to zero
cyp2c9_30x_insert[is.na(cyp2c9_30x_insert)] <- 0

## Combine insert size data from all coverages ##
cyp2c9_insert <- rbind(cyp2c9_5x_insert, cyp2c9_15x_insert, cyp2c9_30x_insert)

# Add unique identifier for each allele x coverage combination
cyp2c9_insert <- cyp2c9_insert %>% mutate("ID" = paste0(Allele, "_", Coverage))

# Make coverage a factor with ordered levels
cyp2c9_insert$Coverage <- factor(cyp2c9_insert$Coverage, levels = c("5x", "15x", "30x"))

## Plot insert size distribution for each allele ##

ggplot(data = cyp2c9_insert, aes(x = Insert_Size, y = Percent)) + 
  geom_line(aes(group = ID, color = Coverage), linewidth = 1, alpha = 0.7) + 
  theme_bw() +
  labs(x = "Insert Size (bp)", y = "Frequency") + xlim(0, 2000) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

## Take average over all simulated alleles (excluding CYP2C9*1)

cyp2c9_insert_grouped <- cyp2c9_insert %>% filter(Allele != "CYP2C9*1") %>%
  group_by(Insert_Size) %>% 
  summarize("Percent" = mean(Percent), .groups = "drop")

cyp2c9_insert_grouped$Data <- "Synthetic"
cyp2c9_insert_grouped$Group <- "Variants"

# Add CYP2C9*1 results

cyp2c9_1_insert <- cyp2c9_insert[cyp2c9_insert$Allele == "CYP2C9*1",] %>%
  group_by(Insert_Size) %>%
  summarize("Percent" = mean(Percent)) %>%
  mutate(Group = "Wild type", Data = "Synthetic")
cyp2c9_insert_grouped <- rbind(cyp2c9_insert_grouped, cyp2c9_1_insert)

# Add real NovaSeq data results
cyp2c9_real_insert_files <- list.files("evaluation_output", 
                                       pattern = "*insert_sizes.txt", full.names=TRUE)

for (sample in 1:length(cyp2c9_real_insert_files)) {
  
  insert_file <- cyp2c9_real_insert_files[sample]
  insert_data <- read.table(insert_file, sep = "\t", col.names = c("Insert_Size", "Count"))
  
  # Obtain counts and percentages of each insert size
  IS_count <- insert_data$Count
  insert_data$Percent <- IS_count / sum(IS_count)
  insert_sizes <- insert_data$Insert_Size
  
  # Assign label
  label <- rep("Real", times = length(insert_sizes_all))
  
  # Initiate data frame 
  if (sample == 1) {
    cyp2c9_real_insert_all <- merge(IS_data_join, insert_data, by = "Insert_Size", all.x = TRUE)
    cyp2c9_real_insert_all$Data <- label
    cyp2c9_real_insert_all$Group <- label

  }
  else {
    data_temp <- merge(IS_data_join, insert_data, by = "Insert_Size", all.x = TRUE)
    data_temp$Data <- label
    data_temp$Group <- label
    
    cyp2c9_real_insert_all <- rbind(cyp2c9_real_insert_all, data_temp)
  }
}

# Set all NA values to zero
cyp2c9_real_insert_all[is.na(cyp2c9_real_insert_all)] <- 0

# Take average of frequencies over all samples 

cyp2c9_real_insert <- cyp2c9_real_insert_all %>% group_by(Insert_Size) %>% 
  summarize("Percent" = mean(Percent), .groups = "drop") %>%
  mutate(Data = "Real", Group = "Real")

# Add real data results to simulated data results

cyp2c9_insert_grouped <- rbind(cyp2c9_insert_grouped, cyp2c9_real_insert)

## Plot insert size distribution for wild type and variants (average) ##

ggplot(data = cyp2c9_insert_grouped, aes(x = Insert_Size, y = Percent)) + 
  geom_line(aes(group = Group, color = Group, linetype = Data), 
            linewidth = 1.5, alpha = 0.7) + 
  theme_bw() + labs(x = "Insert Size (bp)", y = "Frequency") + xlim(0, 2000) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

## Take average over all simulated alleles, including CYP2C9*1

cyp2c9_insert_grouped2 <- cyp2c9_insert %>% group_by(Insert_Size) %>% 
  summarize("Percent" = mean(Percent), .groups = "drop")

cyp2c9_insert_grouped2$Data <- "Synthetic"

# Add real NovaSeq data results
cyp2c9_insert_grouped2 <- rbind(cyp2c9_insert_grouped2, cyp2c9_real_insert %>%
                                  dplyr::select(!Group))

## Plot insert size distribution for simulated data and real data ##

p_insert <- ggplot(data = cyp2c9_insert_grouped2, aes(x = Insert_Size, y = Percent)) + 
  geom_line(aes(group = Data, color = Data), 
            linewidth = 1.5, alpha = 0.7) + 
  theme_bw() + labs(x = "Insert Size (bp)", y = "Frequency") + xlim(0, 2000) +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

# Make vectors expanding GC content for real and simulated data
cyp2c9_insert_total_counts <- rbind(cyp2c9_insert %>% dplyr::select(c(Insert_Size, Count, Data)),
                                cyp2c9_real_insert_all %>% dplyr::select(c(Insert_Size, Count, Data)))

cyp2c9_insert_total_counts <- cyp2c9_insert_total_counts %>% group_by(Data, Insert_Size) %>%
  summarize("Total-Count" = sum(Count), .groups = "drop")

insert_sizes_synthetic_reads <- c()
insert_sizes_real_reads <- c()

for (i in insert_sizes_all) {
  
  # Get real data count
  count_real <- (cyp2c9_insert_total_counts %>% 
                   filter(Data == "Real", Insert_Size == i))$`Total-Count`
  
  # Get simulated data count
  count_synthetic <- (cyp2c9_insert_total_counts %>% 
                        filter(Data == "Synthetic", Insert_Size == i))$`Total-Count`
  
  # Add to vectors
  insert_sizes_real_reads <- c(insert_sizes_real_reads, rep(i, times = count_real))
  insert_sizes_synthetic_reads <- c(insert_sizes_synthetic_reads, rep(i, times = count_synthetic))
  
}

# Compute means directly
mean(insert_sizes_real_reads)
mean(insert_sizes_synthetic_reads)

# Compute variances directly
insert_var_real <- var(insert_sizes_real_reads)
insert_var_synthetic <- var(insert_sizes_synthetic_reads)

sqrt(insert_var_real)
sqrt(insert_var_synthetic)


##----------------------------------------------------------------------------##

### COVERAGE ###

# -> without comparison to real data!

##----------------------------------------------------------------------------##

## Save:

ggarrange(p_qscores, p_error, nrow = 1, ncol = 2, common.legend = TRUE,
          legend = "bottom")

ggarrange(p_GC, p_insert, nrow = 1, ncol = 2, common.legend = TRUE,
          legend = "bottom")

ggarrange(p_cyp2c9_qscores_read1, p_cyp2c9_qscores_read2, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

ggarrange(p_cyp2c9_qscores_grouped_read1, p_cyp2c9_qscores_grouped_read2, 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

ggarrange(p_sub, p_indel, nrow = 1, ncol = 2, common.legend = TRUE,
          legend = "bottom")

ggarrange(p_cyp2c9_err_grouped_read1, p_cyp2c9_err_grouped_read2, 
          nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")

ggarrange(p_cyp2c9_err_read1, p_cyp2c9_err_read2, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")

