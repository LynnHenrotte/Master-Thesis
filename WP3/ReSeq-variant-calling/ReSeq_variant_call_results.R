###-------------------------------###
### Master Thesis - Lynn Henrotte ###
###-------------------------------############################################
### Visualizations for variant calling quality scores of ReSeq simulations ###
###------------------------------------------------------------------------###

## Packages ##

library(tidyverse)
library(ggplot2)
library(ggpubr)

## Preliminaries ##

path_to_folder = paste0("C:/Users/lynnh/OneDrive/Bureaublad",
                        "/2nd Master Stat/Master Thesis/WP3/ReSeq variant call results")
setwd(path_to_folder)

##-------------------------------------------------------------------------##
##---- CYP2C9 -------------------------------------------------------------##
##-------------------------------------------------------------------------##

## Obtain data ##

## Variant calling stats
sim_stats_5x <- read.delim(file = "ReSeq_CYP2C9_5x_sim_stats.tsv", sep = "\t")
sim_stats_15x <- read.delim(file = "ReSeq_CYP2C9_15x_sim_stats.tsv", sep = "\t")
sim_stats_30x <- read.delim(file = "ReSeq_CYP2C9_30x_sim_stats.tsv", sep = "\t")

# Change colnames
colnames(sim_stats_5x) <- tolower(colnames(sim_stats_5x))
colnames(sim_stats_5x)[7] <- "systematic_errors" # = percentage of non-allelic called variants that were found as systematic errors
colnames(sim_stats_15x) <- tolower(colnames(sim_stats_15x))
colnames(sim_stats_15x)[7] <- "systematic_errors"
colnames(sim_stats_30x) <- tolower(colnames(sim_stats_30x))
colnames(sim_stats_30x)[7] <- "systematic_errors"

# Adjust number of errors called (include all non-allelic variants)
sim_stats_5x$error_called <- sim_stats_5x$total_called - sim_stats_5x$x.allelic_called
sim_stats_15x$error_called <- sim_stats_15x$total_called - sim_stats_15x$x.allelic_called
sim_stats_30x$error_called <- sim_stats_30x$total_called - sim_stats_30x$x.allelic_called

# Add all variant calling results together
sim_stats_5x$coverage <- "5x"
sim_stats_15x$coverage <- "15x"
sim_stats_30x$coverage <- "30x"

sim_stats_cyp2c9 <- rbind(sim_stats_5x, sim_stats_15x, sim_stats_30x)
sim_stats_cyp2c9$coverage <- factor(sim_stats_cyp2c9$coverage, levels = c("5x", "15x", "30x"))

# Set NA to zero
sim_stats_cyp2c9[is.na(sim_stats_cyp2c9)] <- 0

## Quality scores
allele_qscores <- read.delim(file = "Allelic_variant_qscores.tsv", sep = "\t")
error_qscores <- read.delim(file = "error_variant_qscores.tsv", sep = "\t")

# Order coverage levels
allele_qscores$Coverage <- factor(allele_qscores$Coverage, levels = c("5x", "15x", "30x"))
error_qscores$Coverage <- factor(error_qscores$Coverage, levels = c("5x", "15x", "30x"))

# Add allelic and error quality scores together
allele_qscores$Type <- "Allelic"
error_qscores$Type <- "Error"
all_cyp2c9_qscores <- rbind(allele_qscores, error_qscores)

## Compute summary statistics ##

# Number and percentage of missed allelic variants, by coverage
sim_stats_cyp2c9 %>% group_by(coverage) %>% 
  summarize("num_allelic_missed" = sum(x.allelic_variants - x.allelic_called),
            "percent_allelic_missed" = 100*sum(x.allelic_variants - x.allelic_called)/sum(x.allelic_variants))

# Average number of variants called, by coverage
sim_stats_cyp2c9 %>% group_by(coverage) %>% summarize("mean_total_called" = mean(total_called))

# Average quality score of allelic and non-allelic variants, by coverage
all_cyp2c9_qscores %>% group_by(Type, Coverage) %>% 
  summarize("mean_qscore" = mean(QScore), 
            "sd_qscore" = sd(QScore), .groups = "drop")

# Number of alleles for which an erroneous variant call had higher quality than a correct variant call
sim_stats_cyp2c9 %>% group_by(coverage) %>% 
  summarize("error_Q_gt_allele_Q" = sum(max_qual_errors > min_qual_allelic)) %>%
  mutate("error_Q_gt_allele_Q_percent" = error_Q_gt_allele_Q/85*100)

## Make plots ##

# Density of quality scores of allelic and error-induced variants, at 30x coverage
p_qual_30x <- ggplot(data = all_cyp2c9_qscores %>% filter(Coverage == "30x"), 
        aes(x = QScore, fill = Type)) + theme_bw() + xlim(0, 1780) + ylim(0,0.03) +
  geom_density(adjust = 1, alpha = 0.7) + ggtitle("30x coverage") +
  labs(x = "Variant calling quality score", fill = "Variant Type") +
  scale_fill_manual(labels = c("Allelic", "Error-induced"),
                    values = c("darkolivegreen4", "violetred4")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

# Density of quality scores of allelic and error-induced variants, at 15x coverage
p_qual_15x <- ggplot(data = all_cyp2c9_qscores %>% filter(Coverage == "15x"), 
                     aes(x = QScore, fill = Type)) + theme_bw() + xlim(0, 1780) + ylim(0,0.03) +
  geom_density(adjust = 1, alpha = 0.7) + ggtitle("15x coverage") +
  labs(x = "Variant calling quality score", fill = "Variant Type") +
  scale_fill_manual(labels = c("Allelic", "Error-induced"),
                    values = c("darkolivegreen4", "violetred4")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

# Density of quality scores of allelic and error-induced variants, at 5x coverage
p_qual_5x <- ggplot(data = all_cyp2c9_qscores %>% filter(Coverage == "5x"), 
                     aes(x = QScore, fill = Type)) + theme_bw() + xlim(0, 1780) + ylim(0,0.03) +
  geom_density(adjust = 1, alpha = 0.7) + ggtitle("5x coverage") +
  labs(x = "Variant calling quality score", fill = "Variant Type") +
  scale_fill_manual(labels = c("Allelic", "Error-induced"),
                    values = c("darkolivegreen4", "violetred4")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 25, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

# Combine densities
ggarrange(p_qual_5x, p_qual_15x, p_qual_30x, common.legend = TRUE,
          legend = "bottom", nrow = 1, ncol = 3)

# Density of qscores of allelic variants and maximal qscores of error-induced variants, at 30x coverage
data_max_error_30x <- sim_stats_cyp2c9 %>% filter(coverage == "30x") %>%
  dplyr::select(max_qual_errors) %>% 
  mutate("Coverage" = "30x", "QScore" = max_qual_errors, "Type" = "Error") %>%
  select(!c(max_qual_errors))

data_max_error_30x <- rbind(data_max_error_30x, 
                            all_cyp2c9_qscores %>% filter(Coverage == "30x",
                                                          Type == "Allelic"))
data_max_error_30x$Type <- factor(data_max_error_30x$Type, levels = c("Error", "Allelic"))

p_max_qual_error_30x <- ggplot(data = data_max_error_30x, 
                     aes(x = QScore, fill = Type)) + theme_bw() + xlim(0, 2000) +
  geom_density(adjust = 0.8, alpha = 0.7) + ggtitle("30x coverage") + ylim(0, 0.02) +
  labs(x = "Variant calling quality score", fill = "Variant Type") +
  scale_fill_manual(labels = c("Allelic", "Maximal error-induced"),
                    values = c("darkolivegreen3", "violetred3")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))  

p_max_qual_error_30x

# Density of qscores of allelic variants and maximal qscores of error-induced variants, at 15x coverage
data_max_error_15x <- sim_stats_cyp2c9 %>% filter(coverage == "15x") %>%
  dplyr::select(max_qual_errors) %>% 
  mutate("Coverage" = "15x", "QScore" = max_qual_errors, "Type" = "Error") %>%
  select(!c(max_qual_errors))

data_max_error_15x <- rbind(data_max_error_15x, 
                            all_cyp2c9_qscores %>% filter(Coverage == "15x",
                                                          Type == "Allelic"))
data_max_error_15x$Type <- factor(data_max_error_15x$Type, levels = c("Error", "Allelic"))

p_max_qual_error_15x <- ggplot(data = data_max_error_15x, 
                               aes(x = QScore, fill = Type)) + theme_bw() + xlim(0, 2000) +
  geom_density(adjust = 0.8, alpha = 0.7) + ggtitle("15x coverage") + ylim(0, 0.02) +
  labs(x = "Variant calling quality score", fill = "Variant Type") +
  scale_fill_manual(labels = c("Allelic", "Maximal error-induced"),
                    values = c("darkolivegreen3", "violetred3")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))  

p_max_qual_error_15x

# Density of qscores of allelic variants and maximal qscores of error-induced variants, at 5x coverage
data_max_error_5x <- sim_stats_cyp2c9 %>% filter(coverage == "5x") %>%
  dplyr::select(max_qual_errors) %>% 
  mutate("Coverage" = "5x", "QScore" = max_qual_errors, "Type" = "Error") %>%
  select(!c(max_qual_errors))

data_max_error_5x <- rbind(data_max_error_5x, 
                            all_cyp2c9_qscores %>% filter(Coverage == "5x",
                                                          Type == "Allelic"))
data_max_error_5x$Type <- factor(data_max_error_5x$Type, levels = c("Error", "Allelic"))

p_max_qual_error_5x <- ggplot(data = data_max_error_5x, 
                               aes(x = QScore, fill = Type)) + theme_bw() + xlim(0, 2000) +
  geom_density(adjust = 0.8, alpha = 0.7) + ggtitle("5x coverage") + ylim(0, 0.02) +
  labs(x = "Variant calling quality score", fill = "Variant Type") +
  scale_fill_manual(labels = c("Allelic", "Maximal error-induced"),
                    values = c("darkolivegreen3", "violetred3")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))  

p_max_qual_error_5x

# Combine densities
ggarrange(p_max_qual_error_5x, p_max_qual_error_15x, p_max_qual_error_30x,
          common.legend = TRUE, legend = "bottom", nrow = 1, ncol = 3)

# Histogram of qscores of allelic variants and maximal qscores of error-induced variants, at 30x coverage
hist_max_qual_error_30x <- ggplot(data = data_max_error_30x, 
                               aes(x = QScore, fill = Type)) + theme_bw() + 
  geom_histogram(binwidth = 35, alpha = 0.7) + ggtitle("30x coverage") + coord_cartesian(xlim=c(0, 1800), ylim=c(0, 60)) +
  labs(x = "Variant calling quality score", y = "Frequency", fill = "Variant Type") +
  scale_fill_manual(labels = c("Maximal error-induced", "Allelic"),
                    values = c("violetred3", "darkolivegreen3")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))  

hist_max_qual_error_30x

# Histogram of qscores of allelic variants and maximal qscores of error-induced variants, at 15x coverage
hist_max_qual_error_15x <- ggplot(data = data_max_error_15x, 
                               aes(x = QScore, fill = Type)) + theme_bw() + 
  geom_histogram(binwidth = 35, alpha = 0.7) + ggtitle("15x coverage") + coord_cartesian(xlim=c(0, 1800), ylim=c(0, 60)) +
  labs(x = "", y = "Frequency", fill = "Variant Type") +
  scale_fill_manual(labels = c("Maximal error-induced", "Allelic"),
                    values = c("violetred3", "darkolivegreen3")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))  

hist_max_qual_error_15x

# Density of qscores of allelic variants and maximal qscores of error-induced variants, at 5x coverage
hist_max_qual_error_5x <- ggplot(data = data_max_error_5x, 
                              aes(x = QScore, fill = Type)) + theme_bw() +
  geom_histogram(binwidth = 35, alpha = 0.7) + ggtitle("5x coverage") + coord_cartesian(xlim=c(0, 1800), ylim=c(0, 60)) +
  labs(x = "", y = "Frequency", fill = "Variant Type") +
  scale_fill_manual(labels = c("Maximal error-induced", "Allelic"),
                    values = c("violetred3", "darkolivegreen3")) +
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 20, face = "bold.italic"),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))  

hist_max_qual_error_5x

# Combine histograms
ggarrange(hist_max_qual_error_5x, hist_max_qual_error_15x, hist_max_qual_error_30x,
          common.legend = TRUE, legend = "bottom", nrow = 3, ncol = 1)


# Histogram comparing mean quality scores of allelic variants and error-induced variants
qual_scores <- c(data$MEAN_QUAL_ALLELIC, data$MEAN_QUAL_ERRORS)
groups <- rep(c("Allelic", "Error-induced"), each = 84)
data_hist <- data.frame(qual_scores, groups)

ggplot(data = data_hist, aes(x = qual_scores, fill = groups)) + 
  theme_bw() + geom_histogram(bins = 21) + labs(x = "Average quality score",
                                       y = "Frequency", fill = "Variant") +
  ggtitle("Distribution of average quality scores for allelic and error-induced called variants") +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(breaks = seq(0, 2000, by = 250))


# Histogram comparing the maximal quality score of error-induced variants and the 
# minimal quality score of the allelic variants
qual_scores2 <- c(data$MIN_QUAL_ALLELIC, data$MAX_QUAL_ERRORS)
groups2 <- rep(c("Minimal allelic", "Maximal error-induced"), each = 84)
data_hist2 <- data.frame(qual_scores2, groups2)

ggplot(data = data_hist2, aes(x = qual_scores2, fill = groups2)) + 
  theme_bw() + geom_density(alpha = 0.7) + labs(x = "Quality score",
                                                y = "Frequency", fill = "Quality score") +
  ggtitle("Distribution of minimal/maximal quality scores of allelic/error-induced called variants") +
  theme(text = element_text(size = 20)) + xlim(0,2000) +
    scale_fill_manual(values = c("lightgreen", "lightblue"))


# Histogram of differences in minimal quality score of allelic variants and maximal
# quality score of error-induced variants
data$SMALLEST_QUAL_DIFF <- data$MIN_QUAL_ALLELIC - data$MAX_QUAL_ERR
ggplot(data, aes(x = SMALLEST_QUAL_DIFF)) + theme_bw() + 
  geom_histogram(col = "coral", fill = "coral", alpha = 0.8) +
  geom_vline(xintercept = 0, col = "red2", lty = 2, size = 1) +
  ggtitle("Histogram of smallest quality differences between allelic and error-induced variants") +
  theme(text = element_text(size = 20)) + 
  labs(x = "Smallest quality score difference", y = "Frequency") 
  








