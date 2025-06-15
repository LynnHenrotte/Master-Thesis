### Master Thesis: synthetic read simulation for PGx ###
### Lynn Henrotte 2024-2025 ###

### Diplotype calling: final plots ###

### Packages ###

library("tidyverse")
library("ggplot2")
library("ggpubr")
library("RColorBrewer")
library("VennDiagram")

####------------####
####-- CYP2C9 --####
####------------####

####-- PyPGx --####

#### Barplot indicating how correct the PyPGx diplotype calls were, split up by coverage ####

correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), times = 3), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
coverage <- factor(rep(c("10x", "30x", "60x"), each = 3))
number_correct <- c(5, 2015, 536, 30, 2050, 476, 3, 2122, 413)
percent_correct <- number_correct / 2556
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_cyp2c9_pypgx_coverage <- data.frame(correct, number_correct, percent_correct, 
                                        coverage, percent_correct_label)

bar_pypgx_cyp2c9_coverage <- ggplot(data = data_cyp2c9_pypgx_coverage, 
                                     aes(x = coverage, y = percent_correct,
                                         color = correct, fill = correct,
                                         label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct PyPGx CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom")

bar_pypgx_cyp2c9_coverage2 <- ggplot(data = data_cyp2c9_pypgx_coverage, 
                                      aes(x = correct, y = percent_correct,
                                          color = coverage, fill = coverage,
                                          label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.07) +
  labs(x = "", y = "Frequency", color = "Coverage", fill = "Coverage") +
  #ggtitle("Frequencies of wrong, partial and fully correct PyPGx CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 7,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("pink3", "skyblue2", "seagreen3")) +
  #scale_color_manual(values = c("pink3", "skyblue2", "seagreen3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") 

bar_pypgx_cyp2c9_coverage
bar_pypgx_cyp2c9_coverage2

#### Barplot indicating how often each CYP2C9 star allele was correctly called by PyPGx at 60x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
pypgx_cyp2c9_alleles_60x <- read.table("PyPGx_correct_CYP2C9_alleles_60x_freq.txt",
                                        sep = "\t", header = TRUE)

# Add zero frequency for *25
pypgx_cyp2c9_alleles_60x <- rbind(pypgx_cyp2c9_alleles_60x, c(0, 25))
pypgx_cyp2c9_alleles_60x <- pypgx_cyp2c9_alleles_60x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", pypgx_cyp2c9_alleles_60x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
pypgx_cyp2c9_alleles_60x$Allele_labels <- allele_labels

# Make color factor
pypgx_allele_percent_60x_cyp2c9 <- pypgx_cyp2c9_alleles_60x$Frequency/72
pypgx_cyp2c9_alleles_60x$Percent <- pypgx_allele_percent_60x_cyp2c9
pypgx_allele_colors_60x_cyp2c9 <- ifelse(pypgx_allele_percent_60x_cyp2c9 > 0.85, "Good",
                            ifelse (pypgx_allele_percent_60x_cyp2c9 < 0.4, "Bad", "Intermediate"))
pypgx_allele_colors_60x_cyp2c9 <- factor(pypgx_allele_colors_60x_cyp2c9, 
                            levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
pypgx_allele_colors2_60x_cyp2c9 <- ifelse(pypgx_allele_percent_60x_cyp2c9 >= 0.5, "50% or more", "Less than 50%")
pypgx_allele_colors2_60x_cyp2c9 <- factor(pypgx_allele_colors2_60x_cyp2c9, levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_pypgx_cyp2c9_alleles_60x_1 <- ggplot(data = pypgx_cyp2c9_alleles_60x, 
                                                  aes(x = Allele_labels, y = Percent,
                                                      fill = pypgx_allele_colors_60x_cyp2c9,
                                                      color = pypgx_allele_colors_60x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 60x coverage")

# Create plot with second color factor

bar_pypgx_cyp2c9_alleles_60x_2 <- ggplot(data = pypgx_cyp2c9_alleles_60x, 
                                                  aes(x = Allele_labels, y = Percent,
                                                      fill = pypgx_allele_colors2_60x_cyp2c9,
                                                      color = pypgx_allele_colors2_60x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 60x coverage")

#### Barplot indicating how often each CYP2C9 star allele was correctly called by PyPGx at 30x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
pypgx_cyp2c9_alleles_30x <- read.table("PyPGx_correct_CYP2C9_alleles_30x_freq.txt",
                                       sep = "\t", header = TRUE)

# Add zero frequency for *25
pypgx_cyp2c9_alleles_30x <- rbind(pypgx_cyp2c9_alleles_30x, c(0, 25))
pypgx_cyp2c9_alleles_30x <- pypgx_cyp2c9_alleles_30x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", pypgx_cyp2c9_alleles_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
pypgx_cyp2c9_alleles_30x$Allele_labels <- allele_labels

# Make color factor
pypgx_allele_percent_30x_cyp2c9 <- pypgx_cyp2c9_alleles_30x$Frequency/72
pypgx_cyp2c9_alleles_30x$Percent <- pypgx_allele_percent_30x_cyp2c9
pypgx_allele_colors_30x_cyp2c9 <- ifelse(pypgx_allele_percent_30x_cyp2c9 > 0.85, "Good",
                                   ifelse (pypgx_allele_percent_30x_cyp2c9 < 0.4, "Bad", "Intermediate"))
pypgx_allele_colors_30x_cyp2c9 <- factor(pypgx_allele_colors_30x_cyp2c9, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
pypgx_allele_colors2_30x_cyp2c9 <- ifelse(pypgx_allele_percent_30x_cyp2c9 >= 0.5,
                                          "50% or more", "Less than 50%")
pypgx_allele_colors2_30x_cyp2c9 <- factor(pypgx_allele_colors2_30x_cyp2c9, 
                                          levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_pypgx_cyp2c9_alleles_30x_1 <- ggplot(data = pypgx_cyp2c9_alleles_30x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = pypgx_allele_colors_30x_cyp2c9,
                                             color = pypgx_allele_colors_30x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 30x coverage")

# Create plot with second color factor

bar_pypgx_cyp2c9_alleles_30x_2 <- ggplot(data = pypgx_cyp2c9_alleles_30x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = pypgx_allele_colors2_30x_cyp2c9,
                                             color = pypgx_allele_colors2_30x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 30x coverage")

#### Barplot indicating how often each CYP2C9 star allele was correctly called by PyPGx at 10x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
pypgx_cyp2c9_alleles_10x <- read.table("PyPGx_correct_CYP2C9_alleles_10x_freq.txt",
                                       sep = "\t", header = TRUE)

# Add zero frequency for *25
pypgx_cyp2c9_alleles_10x <- rbind(pypgx_cyp2c9_alleles_10x, c(0, 25))
pypgx_cyp2c9_alleles_10x <- pypgx_cyp2c9_alleles_10x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", pypgx_cyp2c9_alleles_10x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
pypgx_cyp2c9_alleles_10x$Allele_labels <- allele_labels

# Make color factor
pypgx_allele_percent_10x_cyp2c9 <- pypgx_cyp2c9_alleles_10x$Frequency/72
pypgx_cyp2c9_alleles_10x$Percent <- pypgx_allele_percent_10x_cyp2c9
pypgx_allele_colors_10x_cyp2c9 <- ifelse(pypgx_allele_percent_10x_cyp2c9 > 0.85, "Good",
                                   ifelse (pypgx_allele_percent_10x_cyp2c9 < 0.4, "Bad", "Intermediate"))
pypgx_allele_colors_10x_cyp2c9 <- factor(pypgx_allele_colors_10x_cyp2c9, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
pypgx_allele_colors2_10x_cyp2c9 <- ifelse(pypgx_allele_percent_10x_cyp2c9 >= 0.5, 
                                          "50% or more", "Less than 50%")
pypgx_allele_colors2_10x_cyp2c9 <- factor(pypgx_allele_colors2_10x_cyp2c9, 
                                          levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_pypgx_cyp2c9_alleles_10x_1 <- ggplot(data = pypgx_cyp2c9_alleles_10x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = pypgx_allele_colors_10x_cyp2c9,
                                             color = pypgx_allele_colors_10x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 10x coverage")

# Create plot with second color factor

bar_pypgx_cyp2c9_alleles_10x_2 <- ggplot(data = pypgx_cyp2c9_alleles_10x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = pypgx_allele_colors2_10x_cyp2c9,
                                             color = pypgx_allele_colors2_10x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 10x coverage")


## Combine barplots for all coverages

ggarrange(bar_pypgx_cyp2c9_alleles_10x_1,
          bar_pypgx_cyp2c9_alleles_30x_1, 
          bar_pypgx_cyp2c9_alleles_60x_1,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_pypgx_cyp2c9_alleles_10x_2,
          bar_pypgx_cyp2c9_alleles_30x_2, 
          bar_pypgx_cyp2c9_alleles_60x_2,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

####-- Aldy --####

#### Barplot indicating how correct the Aldy diplotype calls were, split up by coverage ####

correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), times = 3), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
coverage <- factor(rep(c("10x", "30x", "60x"), each = 3))
number_correct <- c(9, 309, 3337, 2, 33, 3620, 1, 9, 3645)
percent_correct <- number_correct / 3655
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_aldy_cyp2c9_coverage <- data.frame(correct, number_correct, percent_correct, 
                                       coverage, percent_correct_label)

bar_aldy_cyp2c9_coverage <- ggplot(data = data_aldy_cyp2c9_coverage, 
                                    aes(x = coverage, y = percent_correct,
                                        color = correct, fill = correct,
                                        label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct Aldy CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom")


bar_aldy_cyp2c9_coverage2 <- ggplot(data = data_aldy_cyp2c9_coverage, 
                                     aes(x = correct, y = percent_correct,
                                         color = coverage, fill = coverage,
                                         label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.07) +
  labs(x = "", y = "Frequency", color = "Coverage", fill = "Coverage") +
  #ggtitle("Frequencies of wrong, partial and fully correct Aldy CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.1, size = 7,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("pink3", "skyblue2", "seagreen3")) +
  #scale_color_manual(values = c("pink3", "skyblue2", "seagreen3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom")

bar_aldy_cyp2c9_coverage
bar_aldy_cyp2c9_coverage2

#### Barplot indicating how often each CYP2C9 star allele was correctly called by Aldy at 60x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
aldy_cyp2c9_alleles_60x <- read.table("Aldy_correct_CYP2C9_alleles_60x_freq.txt",
                                       sep = "\t", header = TRUE)

# Make labels
allele_labels <- paste0("*", aldy_cyp2c9_alleles_60x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
aldy_cyp2c9_alleles_60x$Allele_labels <- allele_labels

# Make color factor
aldy_allele_percent_60x_cyp2c9 <- aldy_cyp2c9_alleles_60x$Frequency/86
aldy_cyp2c9_alleles_60x$Percent <- aldy_allele_percent_60x_cyp2c9
aldy_allele_colors_60x_cyp2c9 <- ifelse(aldy_allele_percent_60x_cyp2c9 > 0.85, "Good",
                                   ifelse (aldy_allele_percent_60x_cyp2c9 < 0.4, "Bad", "Intermediate"))
aldy_allele_colors_60x_cyp2c9 <- factor(aldy_allele_colors_60x_cyp2c9, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
aldy_allele_colors2_60x_cyp2c9 <- ifelse(aldy_allele_percent_60x_cyp2c9 >= 0.5, 
                                         "50% or more", "Less than 50%")
aldy_allele_colors2_60x_cyp2c9 <- factor(aldy_allele_colors2_60x_cyp2c9, 
                                         levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_aldy_cyp2c9_alleles_60x_1 <- ggplot(data = aldy_cyp2c9_alleles_60x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = aldy_allele_colors_60x_cyp2c9,
                                             color = aldy_allele_colors_60x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by Aldy at 60x coverage")

# Create plot with second color factor

bar_aldy_cyp2c9_alleles_60x_2 <- ggplot(data = aldy_cyp2c9_alleles_60x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = aldy_allele_colors2_60x_cyp2c9,
                                             color = aldy_allele_colors2_60x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by Aldy at 60x coverage")

#### Barplot indicating how often each CYP2C9 star allele was correctly called by Aldy at 30x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
aldy_cyp2c9_alleles_30x <- read.table("Aldy_correct_CYP2C9_alleles_30x_freq.txt",
                                       sep = "\t", header = TRUE)

# Make labels
allele_labels <- paste0("*", aldy_cyp2c9_alleles_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
aldy_cyp2c9_alleles_30x$Allele_labels <- allele_labels

# Make color factor
aldy_allele_percent_30x_cyp2c9 <- aldy_cyp2c9_alleles_30x$Frequency/86
aldy_cyp2c9_alleles_30x$Percent <- aldy_allele_percent_30x_cyp2c9
aldy_allele_colors_30x_cyp2c9 <- ifelse(aldy_allele_percent_30x_cyp2c9 > 0.85, "Good",
                                   ifelse (aldy_allele_percent_30x_cyp2c9 < 0.4, "Bad", "Intermediate"))
aldy_allele_colors_30x_cyp2c9 <- factor(aldy_allele_colors_30x_cyp2c9, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
aldy_allele_colors2_30x_cyp2c9 <- ifelse(aldy_allele_percent_30x_cyp2c9 >= 0.5, "50% or more", "Less than 50%")
aldy_allele_colors2_30x_cyp2c9 <- factor(aldy_allele_colors2_30x_cyp2c9, levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_aldy_cyp2c9_alleles_30x_1 <- ggplot(data = aldy_cyp2c9_alleles_30x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = aldy_allele_colors_30x_cyp2c9,
                                             color = aldy_allele_colors_30x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by Aldy at 30x coverage")

# Create plot with second color factor

bar_aldy_cyp2c9_alleles_30x_2 <- ggplot(data = aldy_cyp2c9_alleles_30x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = aldy_allele_colors2_30x_cyp2c9,
                                             color = aldy_allele_colors2_30x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by Aldy at 30x coverage")

#### Barplot indicating how often each CYP2C9 star allele was correctly called by Aldy at 10x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
aldy_cyp2c9_alleles_10x <- read.table("Aldy_correct_CYP2C9_alleles_10x_freq.txt",
                                       sep = "\t", header = TRUE)

# Make labels
allele_labels <- paste0("*", aldy_cyp2c9_alleles_10x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
aldy_cyp2c9_alleles_10x$Allele_labels <- allele_labels

# Make color factor
aldy_allele_percent_10x_cyp2c9 <- aldy_cyp2c9_alleles_10x$Frequency/86
aldy_cyp2c9_alleles_10x$Percent <- aldy_allele_percent_10x_cyp2c9
aldy_allele_colors_10x_cyp2c9 <- ifelse(aldy_allele_percent_10x_cyp2c9 > 0.85, "Good",
                                   ifelse (aldy_allele_percent_10x_cyp2c9 < 0.4, "Bad", "Intermediate"))
aldy_allele_colors_10x_cyp2c9 <- factor(aldy_allele_colors_10x_cyp2c9, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
aldy_allele_colors2_10x_cyp2c9 <- ifelse(aldy_allele_percent_10x_cyp2c9 >= 0.5, "50% or more", "Less than 50%")
aldy_allele_colors2_10x_cyp2c9 <- factor(aldy_allele_colors2_10x_cyp2c9, levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_aldy_cyp2c9_alleles_10x_1 <- ggplot(data = aldy_cyp2c9_alleles_10x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = aldy_allele_colors_10x_cyp2c9,
                                             color = aldy_allele_colors_10x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by Aldy at 10x coverage")

# Create plot with second color factor

bar_aldy_cyp2c9_alleles_10x_2 <- ggplot(data = aldy_cyp2c9_alleles_10x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = aldy_allele_colors2_10x_cyp2c9,
                                             color = aldy_allele_colors2_10x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by Aldy at 10x coverage")

## Combine barplots for all coverages

ggarrange(bar_aldy_cyp2c9_alleles_10x_1,
          bar_aldy_cyp2c9_alleles_30x_1, 
          bar_aldy_cyp2c9_alleles_60x_1,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_aldy_cyp2c9_alleles_10x_2,
          bar_aldy_cyp2c9_alleles_30x_2, 
          bar_aldy_cyp2c9_alleles_60x_2,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

####-- StellarPGx --####

#### Barplot indicating how correct the StellarPGx diplotype calls were, split up by coverage ####

correct <- factor(rep(c("Failed", "Wrong", "Partially Correct", "Correct"), times = 3), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
coverage <- factor(rep(c("10x", "30x", "60x"), each = 4))
number_correct <- c(5, 148, 1288, 2214, 6, 38, 617, 2994, 3, 75, 120, 3457)
percent_correct <- number_correct / 3655
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_stellarpgx_cyp2c9_coverage <- data.frame(correct, number_correct, percent_correct, 
                                             coverage, percent_correct_label)

bar_stellarpgx_cyp2c9_coverage <- ggplot(data = data_stellarpgx_cyp2c9_coverage, 
                                          aes(x = coverage, y = percent_correct,
                                              color = correct, fill = correct,
                                              label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct StellarPGx CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom")

bar_stellarpgx_cyp2c9_coverage2 <- ggplot(data = data_stellarpgx_cyp2c9_coverage, 
                                           aes(x = correct, y = percent_correct,
                                               color = coverage, fill = coverage,
                                               label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.07) +
  labs(x = "", y = "Frequency", color = "Coverage", fill = "Coverage") +
  #ggtitle("Frequencies of wrong, partial and fully correct StellarPGx CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 7,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  #scale_fill_manual(values = c("pink3", "skyblue2", "seagreen3")) +
  #scale_color_manual(values = c("pink3", "skyblue2", "seagreen3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom")

bar_stellarpgx_cyp2c9_coverage
bar_stellarpgx_cyp2c9_coverage2

### CONCLUSIONS BASED ON THIS GRAPH: AS COVERAGE DECREASES, THE NUMBER OF CORRECT DIPLOTYPE CALLS BY STELLARPGX
### DECREASES AS WELL. HOWEVER, IT IS PROMISING THAT MOST OF THESE MISSED CALLS END UP AS PARTIALLY CORRECT
### INSTEAD OF COMPLETELY WRONG, THE CATEGORY WHICH DOES NOT SEEM TO CHANGE A LOT WHEN COVERAGE CHANGES.

#### Barplot indicating how often each CYP2C9 star allele was correctly called by StellarPGx at 60x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
stellarpgx_cyp2c9_alleles_60x <- read.table("StellarPGx_correct_CYP2C9_alleles_60x_freq.txt",
                                      sep = "\t", header = TRUE)

# Add 0 frequency for *25
stellarpgx_cyp2c9_alleles_60x <- rbind(stellarpgx_cyp2c9_alleles_60x, c(0, 25))
stellarpgx_cyp2c9_alleles_60x <- stellarpgx_cyp2c9_alleles_60x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", stellarpgx_cyp2c9_alleles_60x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
stellarpgx_cyp2c9_alleles_60x$Allele_labels <- allele_labels

# Make color factor
stellarpgx_allele_percent_60x_cyp2c9 <- stellarpgx_cyp2c9_alleles_60x$Frequency/86
stellarpgx_cyp2c9_alleles_60x$Percent <- stellarpgx_allele_percent_60x_cyp2c9
stellarpgx_allele_colors_60x_cyp2c9 <- ifelse(stellarpgx_allele_percent_60x_cyp2c9 > 0.85, "Good",
                                   ifelse (stellarpgx_allele_percent_60x_cyp2c9 < 0.4, "Bad", "Intermediate"))
stellarpgx_allele_colors_60x_cyp2c9 <- factor(stellarpgx_allele_colors_60x_cyp2c9, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
stellarpgx_allele_colors2_60x_cyp2c9 <- ifelse(stellarpgx_allele_percent_60x_cyp2c9 >= 0.5, 
                                               "50% or more", "Less than 50%")
stellarpgx_allele_colors2_60x_cyp2c9 <- factor(stellarpgx_allele_colors2_60x_cyp2c9, 
                                               levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor
bar_stellarpgx_cyp2c9_alleles_60x_1 <- ggplot(data = stellarpgx_cyp2c9_alleles_60x, 
                                        aes(x = Allele_labels, y = Percent,
                                            fill = stellarpgx_allele_colors_60x_cyp2c9,
                                            color = stellarpgx_allele_colors_60x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by StellarPGx at 60x coverage")

# Create plot with second color factor
bar_stellarpgx_cyp2c9_alleles_60x_2 <- ggplot(data = stellarpgx_cyp2c9_alleles_60x, 
                                        aes(x = Allele_labels, y = Percent,
                                            fill = stellarpgx_allele_colors2_60x_cyp2c9,
                                            color = stellarpgx_allele_colors2_60x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by StellarPGx at 60x coverage")

#### Barplot indicating how often each CYP2C9 star allele was correctly called by StellarPGx at 30x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
stellarpgx_cyp2c9_alleles_30x <- read.table("StellarPGx_correct_CYP2C9_alleles_30x_freq.txt",
                                      sep = "\t", header = TRUE)
# Add 0 frequency for *25
stellarpgx_cyp2c9_alleles_30x <- rbind(stellarpgx_cyp2c9_alleles_30x, c(0, 25))
stellarpgx_cyp2c9_alleles_30x <- stellarpgx_cyp2c9_alleles_30x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", stellarpgx_cyp2c9_alleles_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
stellarpgx_cyp2c9_alleles_30x$Allele_labels <- allele_labels

# Make color factor
stellarpgx_allele_percent_30x_cyp2c9 <- stellarpgx_cyp2c9_alleles_30x$Frequency/86
stellarpgx_cyp2c9_alleles_30x$Percent <- stellarpgx_allele_percent_30x_cyp2c9
stellarpgx_allele_colors_30x_cyp2c9 <- ifelse(stellarpgx_allele_percent_30x_cyp2c9 > 0.85, "Good",
                                   ifelse (stellarpgx_allele_percent_30x_cyp2c9 < 0.4, "Bad", "Intermediate"))
stellarpgx_allele_colors_30x_cyp2c9 <- factor(stellarpgx_allele_colors_30x_cyp2c9, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
stellarpgx_allele_colors2_30x_cyp2c9 <- ifelse(stellarpgx_allele_percent_30x_cyp2c9 >= 0.5, 
                                               "50% or more", "Less than 50%")
stellarpgx_allele_colors2_30x_cyp2c9 <- factor(stellarpgx_allele_colors2_30x_cyp2c9, 
                                               levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor
bar_stellarpgx_cyp2c9_alleles_30x_1 <- ggplot(data = stellarpgx_cyp2c9_alleles_30x, 
                                        aes(x = Allele_labels, y = Percent,
                                            fill = stellarpgx_allele_colors_30x_cyp2c9,
                                            color = stellarpgx_allele_colors_30x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by StellarPGx at 30x coverage")

# Create plot with second color factor
bar_stellarpgx_cyp2c9_alleles_30x_2 <- ggplot(data = stellarpgx_cyp2c9_alleles_30x, 
                                        aes(x = Allele_labels, y = Percent,
                                            fill = stellarpgx_allele_colors2_30x_cyp2c9,
                                            color = stellarpgx_allele_colors2_30x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by StellarPGx at 30x coverage")

#### Barplot indicating how often each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
stellarpgx_cyp2c9_alleles_10x <- read.table("StellarPGx_correct_CYP2C9_alleles_10x_freq.txt",
                                      sep = "\t", header = TRUE)

# Add zero frequency for *25
stellarpgx_cyp2c9_alleles_10x <- rbind(stellarpgx_cyp2c9_alleles_10x, c(0, 25))
stellarpgx_cyp2c9_alleles_10x <- stellarpgx_cyp2c9_alleles_10x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", stellarpgx_cyp2c9_alleles_10x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
stellarpgx_cyp2c9_alleles_10x$Allele_labels <- allele_labels

# Make color factor
stellarpgx_allele_percent_10x_cyp2c9 <- stellarpgx_cyp2c9_alleles_10x$Frequency/86
stellarpgx_cyp2c9_alleles_10x$Percent <- stellarpgx_allele_percent_10x_cyp2c9
stellarpgx_allele_colors_10x_cyp2c9 <- ifelse(stellarpgx_allele_percent_10x_cyp2c9 > 0.85, "Good",
                                   ifelse (stellarpgx_allele_percent_10x_cyp2c9 < 0.4, "Bad", "Intermediate"))
stellarpgx_allele_colors_10x_cyp2c9 <- factor(stellarpgx_allele_colors_10x_cyp2c9, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
stellarpgx_allele_colors2_10x_cyp2c9 <- ifelse(stellarpgx_allele_percent_10x_cyp2c9 >= 0.5, 
                                               "50% or more", "Less than 50%")
stellarpgx_allele_colors2_10x_cyp2c9 <- factor(stellarpgx_allele_colors2_10x_cyp2c9, 
                                               levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor
bar_stellarpgx_cyp2c9_alleles_10x_1 <- ggplot(data = stellarpgx_cyp2c9_alleles_10x, 
                                        aes(x = Allele_labels, y = Percent,
                                            fill = stellarpgx_allele_colors_10x_cyp2c9,
                                            color = stellarpgx_allele_colors_10x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage")

# Create plot with second color factor
bar_stellarpgx_cyp2c9_alleles_10x_2 <- ggplot(data = stellarpgx_cyp2c9_alleles_10x, 
                                        aes(x = Allele_labels, y = Percent,
                                            fill = stellarpgx_allele_colors2_10x_cyp2c9,
                                            color = stellarpgx_allele_colors2_10x_cyp2c9)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage")

## Combine barplots for all coverages
ggarrange(bar_stellarpgx_cyp2c9_alleles_10x_1,
          bar_stellarpgx_cyp2c9_alleles_30x_1, 
          bar_stellarpgx_cyp2c9_alleles_60x_1,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_stellarpgx_cyp2c9_alleles_10x_2,
          bar_stellarpgx_cyp2c9_alleles_30x_2, 
          bar_stellarpgx_cyp2c9_alleles_60x_2,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")


####-- PyPGx, Aldy and StellarPGx combined --####

ggarrange(bar_pypgx_cyp2c9_coverage2,
          bar_aldy_cyp2c9_coverage2,
          bar_stellarpgx_cyp2c9_coverage2,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

###-- Venn diagrams --###

# Get vector of all cyp2c9 diplotypes
cyp2c9_all <- c()
for (allele1 in 1:85) {
  for (allele2 in allele1:85) {
    diplotype <- paste0("*", allele1, "/*", allele2)
    cyp2c9_all <- c(cyp2c9_all, diplotype)
  }
}

## Venn diagram studying overlap between diplotype calls by different callers

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")

for (coverage in c("10x", "30x", "60x")) {
  
  for (correctness in c("0%", "50%", "100%")) {
  
    # Read in data
    pypgx_cyp2c9_alleles <- readLines(paste0("Venn_diagram_data/PyPGx_CYP2C9_", coverage, "_",
                                             correctness, ".txt"))
    aldy_cyp2c9_alleles <- readLines(paste0("Venn_diagram_data/Aldy_CYP2C9_", coverage, "_",
                                            correctness, ".txt"))
    stellarpgx_cyp2c9_alleles <- readLines(paste0("Venn_diagram_data/StellarPGx_CYP2C9_", coverage, "_",
                                                  correctness, ".txt"))
    
    # Ensure uniqueness of diplotypes
    pypgx_cyp2c9_alleles <- unique(pypgx_cyp2c9_alleles)
  
    aldy_cyp2c9_alleles <- unique(aldy_cyp2c9_alleles)
    
    stellarpgx_cyp2c9_alleles <- unique(stellarpgx_cyp2c9_alleles)
    
    # Combine data in list
    venn_data <- list("PyPGx" = pypgx_cyp2c9_alleles, "Aldy" = aldy_cyp2c9_alleles,
                      "StellarPGx" = stellarpgx_cyp2c9_alleles)
    
    # Draw Venn diagram
    title <- ifelse(correctness == "0%", paste0("Wrong CYP2C9 calls"),
                    ifelse(correctness == "100%", paste0("Correct CYP2C9 calls"),
                           paste0("Partially correct CYP2C9 calls")))
    sub_title <- paste0(coverage, " coverage")
    
    file_name <- ifelse(correctness == "0%", paste0("VennDiagrams/venn_CYP2C9_", coverage, "_wrong.png"),
                        ifelse(correctness == "100%", paste0("VennDiagrams/venn_CYP2C9_", coverage, "_correct.png"),
                               paste0("VennDiagrams/venn_CYP2C9_", coverage, "_partially_correct.png")))
    
    if (coverage == "30x" && correctness == "0%") {
      adjust <- list(c(.5,-4.5), c(-.6, -1.3), c(.58,6))
    } 
    else if (coverage == "60x" && correctness == "50%") {
      adjust <- list(c(.58,6), c(.5,-4.5), c(0.2,-3.5))
    }
    else if (coverage == "60x" && correctness == "0%") {
      adjust <- list(c(0,-4.5), c(1, -1), c(0.58, 20))
    }
    else {
      adjust <- list(c(-0.1,-2), c(1, -1), c(0.5,1.8))
    }
    
    venn.diagram(x = venn_data, main = title, main.cex = 2.5, sub = sub_title, sub.cex = 2,
                 filename = file_name, category.names = c('PyPGx', 'Aldy', 'StellarPGx'),
                 cat.default.pos = "text", margin = 0.05, cat.cex = 1.6, cex = 1.5,
                 col=c("magenta4", "maroon2", "skyblue3"), 
                 fill = c("magenta4", "maroon2", "skyblue3"),
                 fontface = 'italic', cat.fontface = "bold", alpha = 0.3, 
                 print.mode = c('raw', 'percent'), imagetype = "png",
                 cat.just = adjust)
  }
}

# cat.just:
# first argument = horizontal placement (negative -> right and positive -> left)
# second argument = vertical placement (negative -> up and positive -> down)

####-------------####
####-- CYP2C19 --####
####-------------####

####-- PyPGx --####

## Barplot indicating how correct the PyPGx diplotype calls were, at 30x coverage ##

correct <- factor(c("Failed", "Wrong", "Partially Correct", "Correct"), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
number_correct <- c(0, 3, 296, 367)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_pypgx_calls_30x_cyp2c19 <- data.frame(correct, number_correct, percent_correct, 
                                           percent_correct_label)

bar_pypgx_correct_30x_cyp2c19 <- ggplot(data = data_pypgx_calls_30x_cyp2c19, 
                                        aes(x = correct, y = percent_correct,
                                            color = correct, fill = correct,
                                            label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx 30x CYP2C19 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.2, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_30x_cyp2c19

#### Barplot indicating how often each CYP2C19 star allele was correctly called by PyPGx at 30x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
pypgx_cyp2c19_alleles_30x <- read.table("PyPGx_correct_CYP2C19_alleles_30x_freq.txt",
                                       sep = "\t", header = TRUE)
pypgx_cyp2c19_alleles_30x <- pypgx_cyp2c19_alleles_30x %>% arrange(Allele)

# Adjust frequencies to exclude sub alleles
new_freq <- pypgx_cyp2c19_alleles_30x$Frequency
new_freq <- ifelse(new_freq > 0, new_freq - 5, new_freq)
pypgx_cyp2c19_alleles_30x$Frequency <- new_freq

# Make labels
allele_labels <- paste0("*", pypgx_cyp2c19_alleles_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
pypgx_cyp2c19_alleles_30x$Allele_labels <- allele_labels

# Make color factor
pypgx_allele_percent_30x_cyp2c19 <- pypgx_cyp2c19_alleles_30x$Frequency/37
pypgx_cyp2c19_alleles_30x$Percent <- pypgx_allele_percent_30x_cyp2c19
pypgx_allele_colors_30x_cyp2c19 <- ifelse(pypgx_allele_percent_30x_cyp2c19 > 0.85, "Good",
                                   ifelse (pypgx_allele_percent_30x_cyp2c19 < 0.4, "Bad", "Intermediate"))
pypgx_allele_colors_30x_cyp2c19 <- factor(pypgx_allele_colors_30x_cyp2c19, 
                                   levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
pypgx_allele_colors2_30x_cyp2c19 <- ifelse(pypgx_allele_percent_30x_cyp2c19 >= 0.5, 
                                           "50% or more", "Less than 50%")
pypgx_allele_colors2_30x_cyp2c19 <- factor(pypgx_allele_colors2_30x_cyp2c19, 
                                           levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_pypgx_cyp2c19_alleles_30x_1 <- ggplot(data = pypgx_cyp2c19_alleles_30x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = pypgx_allele_colors_30x_cyp2c19,
                                             color = pypgx_allele_colors_30x_cyp2c19)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C19 star allele was correctly called by PyPGx at 30x coverage")

# Create plot with second color factor

bar_pypgx_cyp2c19_alleles_30x_2 <- ggplot(data = pypgx_cyp2c19_alleles_30x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = pypgx_allele_colors2_30x_cyp2c19,
                                             color = pypgx_allele_colors2_30x_cyp2c19)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C19 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_cyp2c19_alleles_30x_1
bar_pypgx_cyp2c19_alleles_30x_2

####-- Aldy --####

## Barplot indicating how correct the Aldy diplotype calls were, at 30x coverage ##

correct <- factor(c("Failed", "Wrong", "Partially Correct", "Correct"), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
number_correct <- c(0, 3, 70, 593)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_aldy_calls_30x_cyp2c19 <- data.frame(correct, number_correct, percent_correct, 
                                           percent_correct_label)

bar_aldy_correct_30x_cyp2c19 <- ggplot(data = data_aldy_calls_30x_cyp2c19, 
                                        aes(x = correct, y = percent_correct,
                                            color = correct, fill = correct,
                                            label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy 30x CYP2C19 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.2, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_aldy_correct_30x_cyp2c19

#### Barplot indicating how often each CYP2C19 star allele was correctly called by Aldy at 30x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
aldy_cyp2c19_alleles_30x <- read.table("Aldy_correct_CYP2C19_alleles_30x_freq.txt",
                                        sep = "\t", header = TRUE)
aldy_cyp2c19_alleles_30x <- aldy_cyp2c19_alleles_30x %>% arrange(Allele)

# Adjust frequencies to exclude sub alleles
new_freq <- aldy_cyp2c19_alleles_30x$Frequency
new_freq <- ifelse(new_freq > 0, new_freq - 5, new_freq)
aldy_cyp2c19_alleles_30x$Frequency <- new_freq

# Make labels
allele_labels <- paste0("*", aldy_cyp2c19_alleles_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
aldy_cyp2c19_alleles_30x$Allele_labels <- allele_labels

# Make color factor
aldy_allele_percent_30x_cyp2c19 <- aldy_cyp2c19_alleles_30x$Frequency/37
aldy_cyp2c19_alleles_30x$Percent <- aldy_allele_percent_30x_cyp2c19
aldy_allele_colors_30x_cyp2c19 <- ifelse(aldy_allele_percent_30x_cyp2c19 > 0.85, "Good",
                                    ifelse (aldy_allele_percent_30x_cyp2c19 < 0.4, "Bad", "Intermediate"))
aldy_allele_colors_30x_cyp2c19 <- factor(aldy_allele_colors_30x_cyp2c19, 
                                    levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
aldy_allele_colors2_30x_cyp2c19 <- ifelse(aldy_allele_percent_30x_cyp2c19 >= 0.5, 
                                          "50% or more", "Less than 50%")
aldy_allele_colors2_30x_cyp2c19 <- factor(aldy_allele_colors2_30x_cyp2c19, 
                                          levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_aldy_cyp2c19_alleles_30x_1 <- ggplot(data = aldy_cyp2c19_alleles_30x, 
                                          aes(x = Allele_labels, y = Percent,
                                              fill = aldy_allele_colors_30x_cyp2c19,
                                              color = aldy_allele_colors_30x_cyp2c19)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C19 star allele was correctly called by Aldy at 30x coverage")

# Create plot with second color factor

bar_aldy_cyp2c19_alleles_30x_2 <- ggplot(data = aldy_cyp2c19_alleles_30x, 
                                          aes(x = Allele_labels, y = Percent,
                                              fill = aldy_allele_colors2_30x_cyp2c19,
                                              color = aldy_allele_colors2_30x_cyp2c19)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C19 star allele was correctly called by Aldy at 30x coverage")

bar_aldy_cyp2c19_alleles_30x_1
bar_aldy_cyp2c19_alleles_30x_2

####-- StellarPGx --####

## Barplot indicating how correct the StellarPGx diplotype calls were, at 30x coverage ##

correct <- factor(c("Failed", "Wrong", "Partially Correct", "Correct"), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
number_correct <- c(1, 1, 62, 598)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_stellarpgx_calls_30x_cyp2c19 <- data.frame(correct, number_correct, percent_correct, 
                                           percent_correct_label)

bar_stellarpgx_correct_30x_cyp2c19 <- ggplot(data = data_stellarpgx_calls_30x_cyp2c19, 
                                        aes(x = correct, y = percent_correct,
                                            color = correct, fill = correct,
                                            label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx 30x CYP2C19 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.2, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_stellarpgx_correct_30x_cyp2c19

#### Barplot indicating how often each CYP2C19 star allele was correctly called by StellarPGx at 30x coverage ####

# Read in data
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
stellarpgx_cyp2c19_alleles_30x <- read.table("StellarPGx_correct_CYP2C19_alleles_30x_freq.txt",
                                       sep = "\t", header = TRUE)

# Adjust frequencies for non-consensus calls and add 0 frequency for *37
stellarpgx_cyp2c19_alleles_30x <- rbind(stellarpgx_cyp2c19_alleles_30x, c(0, 37))
stellarpgx_cyp2c19_alleles_30x$Frequency <- ifelse(stellarpgx_cyp2c19_alleles_30x$Allele %in% c(15, 17, 38),
                                                   stellarpgx_cyp2c19_alleles_30x$Frequency + 1,
                                                   stellarpgx_cyp2c19_alleles_30x$Frequency)
stellarpgx_cyp2c19_alleles_30x <- stellarpgx_cyp2c19_alleles_30x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", stellarpgx_cyp2c19_alleles_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
stellarpgx_cyp2c19_alleles_30x$Allele_labels <- allele_labels

# Make color factor
stellarpgx_allele_percent_30x_cyp2c19 <- stellarpgx_cyp2c19_alleles_30x$Frequency/37
stellarpgx_cyp2c19_alleles_30x$Percent <- stellarpgx_allele_percent_30x_cyp2c19
stellarpgx_allele_colors_30x_cyp2c19 <- ifelse(stellarpgx_allele_percent_30x_cyp2c19 > 0.85, "Good",
                                    ifelse (stellarpgx_allele_percent_30x_cyp2c19 < 0.4, "Bad", "Intermediate"))
stellarpgx_allele_colors_30x_cyp2c19 <- factor(stellarpgx_allele_colors_30x_cyp2c19, 
                                    levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
stellarpgx_allele_colors2_30x_cyp2c19 <- ifelse(stellarpgx_allele_percent_30x_cyp2c19 >= 0.5, 
                                                "50% or more", "Less than 50%")
stellarpgx_allele_colors2_30x_cyp2c19 <- factor(stellarpgx_allele_colors2_30x_cyp2c19, 
                                                levels = c("50% or more", "Less than 50%"))

# Create plot with first color factor

bar_stellarpgx_cyp2c19_alleles_30x_1 <- ggplot(data = stellarpgx_cyp2c19_alleles_30x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = stellarpgx_allele_colors_30x_cyp2c19,
                                             color = stellarpgx_allele_colors_30x_cyp2c19)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C19 star allele was correctly called by StellarPGx at 30x coverage")

# Create plot with second color factor

bar_stellarpgx_cyp2c19_alleles_30x_2 <- ggplot(data = stellarpgx_cyp2c19_alleles_30x, 
                                         aes(x = Allele_labels, y = Percent,
                                             fill = stellarpgx_allele_colors2_30x_cyp2c19,
                                             color = stellarpgx_allele_colors2_30x_cyp2c19)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C19 star allele was correctly called by StellarPGx at 30x coverage")

bar_stellarpgx_cyp2c19_alleles_30x_1
bar_stellarpgx_cyp2c19_alleles_30x_2

####-- PyPGx, Aldy, and StellarPGx --####

data_all_calls_30x_cyp2c19 <- rbind(data_pypgx_calls_30x_cyp2c19,
                                    data_aldy_calls_30x_cyp2c19,
                                    data_stellarpgx_calls_30x_cyp2c19)
data_all_calls_30x_cyp2c19$Caller <- factor(rep(c("PyPGx", "Aldy", "StellarPGx"), each = 4),
                                            levels = c("PyPGx", "Aldy", "StellarPGx"))

bar_all_callers_cyp2c19 <- ggplot(data = data_all_calls_30x_cyp2c19, 
                                          aes(x = correct, y = percent_correct,
                                              color = Caller, fill = Caller,
                                              label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.9) + 
  theme_bw() + ylim(0,1) +
  labs(x = "", y = "Frequency", color = "Star allele caller", fill = "Star allele caller") +
  #ggtitle("Frequencies of wrong, partial and fully correct CYP2C19 calls, by different diplotype callers") +
  geom_text(position = position_dodge(width = 0.9), vjust = -.1, size = 7,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text = element_text(size = 25, face = "bold"),
        axis.title = element_text(size = 25, face = "bold.italic"),
        legend.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20, face = "bold.italic"),
        legend.position = "bottom",
        plot.title = element_text(size= 25, hjust = 0.5))

bar_all_callers_cyp2c19

# combine:

ggarrange(bar_pypgx_cyp2c19_alleles_30x_2,
          bar_aldy_cyp2c19_alleles_30x_2,
          bar_stellarpgx_cyp2c19_alleles_30x_2,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

###-- Venn diagrams --###

## Venn diagram studying overlap between CYP2C19 diplotype calls by different callers

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")

for (correctness in c("0%", "50%", "100%")) {
    
  # Read in data
  pypgx_cyp2c19_alleles <- readLines(paste0("Venn_diagram_data/PyPGx_CYP2C19_30x_",
                                           correctness, ".txt"))
  aldy_cyp2c19_alleles <- readLines(paste0("Venn_diagram_data/Aldy_CYP2C19_30x_",
                                          correctness, ".txt"))
  stellarpgx_cyp2c19_alleles <- readLines(paste0("Venn_diagram_data/StellarPGx_CYP2C19_30x_",
                                                correctness, ".txt"))
  
  # Ensure uniqueness of diplotypes
  pypgx_cyp2c19_alleles <- unique(pypgx_cyp2c19_alleles)
  
  aldy_cyp2c19_alleles <- unique(aldy_cyp2c19_alleles)
  
  stellarpgx_cyp2c19_alleles <- unique(stellarpgx_cyp2c19_alleles)
  
  # Combine data in list
  venn_data <- list("PyPGx" = pypgx_cyp2c19_alleles, "Aldy" = aldy_cyp2c19_alleles,
                    "StellarPGx" = stellarpgx_cyp2c19_alleles)
  
  # Draw Venn diagram
  title <- ifelse(correctness == "0%", paste0("Wrong CYP2C19 calls"),
                  ifelse(correctness == "100%", paste0("Correct CYP2C19 calls"),
                         paste0("Partially correct CYP2C19 calls")))
  
  file_name <- ifelse(correctness == "0%", paste0("VennDiagrams/venn_CYP2C19_30x_wrong.png"),
                      ifelse(correctness == "100%", paste0("VennDiagrams/venn_CYP2C19_30x_correct.png"),
                             paste0("VennDiagrams/venn_CYP2C19_30x_partially_correct.png")))
  
  if (correctness == "0%") {
    adjust <- list(c(-0.1,-2), c(1, -1), c(0.5,1.8))
  }
  else if (correctness == "50%") {
    adjust <- list(c(-.5,-7), c(0, -1), c(1,6))
  }
  else {
    adjust <- list(c(1, 5), c(-1, -7), c(.1,6.5))
  }
  
  venn.diagram(x = venn_data, main = title, main.cex = 2.5, sub = "30x coverage", sub.cex = 2,
               filename = file_name, category.names = c('PyPGx', 'Aldy', 'StellarPGx'),
               cat.default.pos = "text", margin = 0.05, cat.cex = 1.6, cex = 1.5,
               col=c("magenta4", "maroon2", "skyblue3"), 
               fill = c("magenta4", "maroon2", "skyblue3"),
               fontface = 'italic', cat.fontface = "bold", alpha = 0.3, 
               print.mode = c('raw', 'percent'), imagetype = "png",
               cat.just = adjust)
}

# cat.just:
# first argument = horizontal placement (negative -> right and positive -> left)
# second argument = vertical placement (negative -> up and positive -> down)

###-------------------------###
###-- OVERALL PERFORMANCE --###
###-------------------------###

# Compute the overall concordance rates for Aldy, PyPGx and StellarPGx, 
# across all alleles of both genes, across all coverages

### PyPGx ###

## Only CYP2C9 ##
pypgx_cyp2c9_60x_correct <- 413
pypgx_cyp2c9_30x_correct <- 476
pypgx_cyp2c9_10x_correct <- 536
(pypgx_cyp2c9_60x_correct + pypgx_cyp2c9_30x_correct + pypgx_cyp2c9_10x_correct)/(3*2556)
# 18.58%

## Only CYP2C19 ##
pypgx_cyp2c19_correct <- 367
pypgx_cyp2c19_correct / 666
# 55.11%

## Both ##
(pypgx_cyp2c9_60x_correct + pypgx_cyp2c9_30x_correct + 
    pypgx_cyp2c9_10x_correct + pypgx_cyp2c19_correct)/(3*2556 + 666)
# 21.50%

### Aldy ###

## Only CYP2C9 ##
aldy_cyp2c9_60x_correct <- 3645
aldy_cyp2c9_30x_correct <- 3620
aldy_cyp2c9_10x_correct <- 3337
(aldy_cyp2c9_60x_correct + aldy_cyp2c9_30x_correct + 
    aldy_cyp2c9_10x_correct)/(3*3655)
# 96.69%

## Only CYP2C19 ##
aldy_cyp2c19_correct <- 593
aldy_cyp2c19_correct/666
# 89.04%

## Both ##
(aldy_cyp2c9_60x_correct + aldy_cyp2c9_30x_correct + 
    aldy_cyp2c9_10x_correct + aldy_cyp2c19_correct)/(3*3655 + 666)
# 96.25%

### StellarPGx ###

## Only CYP2C9 ##
stellarpgx_cyp2c9_60x_correct <- 3457
stellarpgx_cyp2c9_30x_correct <- 2994
stellarpgx_cyp2c9_10x_correct <- 2214
(stellarpgx_cyp2c9_60x_correct + stellarpgx_cyp2c9_30x_correct + 
    stellarpgx_cyp2c9_10x_correct)/(3*3655)
# 79.02%

## Only CYP2C19 ##
stellarpgx_cyp2c19_correct <- 598
stellarpgx_cyp2c19_correct/666
# 89.79%

## Both ##
(stellarpgx_cyp2c9_60x_correct + stellarpgx_cyp2c9_30x_correct + 
    stellarpgx_cyp2c9_10x_correct + stellarpgx_cyp2c19_correct)/(3*3655 + 666)
# 79.64%

###----------------------------###
###-- MOST IMPORTANT ALLELES --###
###----------------------------###

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5/Venn_diagram_data")

# For CYP2C9: *1, *2, *3, *5, *6, *8, *11
correct_cyp2c9_10x_pypgx <- readLines("PyPGx_CYP2C9_10x_100%.txt")
correct_cyp2c9_30x_pypgx <- readLines("PyPGx_CYP2C9_30x_100%.txt")
correct_cyp2c9_60x_pypgx <- readLines("PyPGx_CYP2C9_60x_100%.txt")

correct_cyp2c9_10x_aldy <- readLines("Aldy_CYP2C9_10x_100%.txt")
correct_cyp2c9_30x_aldy <- readLines("Aldy_CYP2C9_30x_100%.txt")
correct_cyp2c9_60x_aldy <- readLines("Aldy_CYP2C9_60x_100%.txt")

correct_cyp2c9_10x_stellarpgx <- readLines("StellarPGx_CYP2C9_10x_100%.txt")
correct_cyp2c9_30x_stellarpgx <- readLines("StellarPGx_CYP2C9_30x_100%.txt")
correct_cyp2c9_60x_stellarpgx <- readLines("StellarPGx_CYP2C9_60x_100%.txt")

cyp2c9_important_alleles <- c(2, 3, 5, 6, 8, 11)
total_cyp2c9_important_diplo <- sum((85-length(cyp2c9_important_alleles)+1):85)

all_cyp2c9_diplo <- c()

for (allele1 in 1:85) {
  for (allele2 in allele1:85) {
    all_cyp2c9_diplo <- c(all_cyp2c9_diplo, paste0("*",allele1,"/*",allele2))
  }
}

# Function that extracts diplotypes containing a certain allele from a vector
get_diplo <- function(input, alleles) {
  
  diplo_out <- c()
  
  for (allele in alleles) {
    
    pattern1 <- paste0("^[*]{1}", allele, "/")
    pattern2 <- paste0("/[*]{1}", allele, "$")
    
    diplo_out1 <- input[grep(pattern1, input, fixed = FALSE, perl = FALSE)]
    diplo_out2 <- input[grep(pattern2, input, fixed = FALSE, perl = FALSE)]
    
    
    diplo_out <- c(diplo_out, diplo_out1, diplo_out2)
  }
  
  return(unique(diplo_out))
  
}

length(unique(get_diplo(all_cyp2c9_diplo, cyp2c9_important_alleles)))

# Overal PyPGx performance on important CYP2C9 alleles 
cyp2c9_10x_pypgx_select <- get_diplo(correct_cyp2c9_10x_pypgx, cyp2c9_important_alleles)
cyp2c9_30x_pypgx_select <- get_diplo(correct_cyp2c9_30x_pypgx, cyp2c9_important_alleles)
cyp2c9_60x_pypgx_select <- get_diplo(correct_cyp2c9_60x_pypgx, cyp2c9_important_alleles)
length(c(cyp2c9_10x_pypgx_select,cyp2c9_30x_pypgx_select,cyp2c9_60x_pypgx_select))/
  (3*total_cyp2c9_important_diplo)
# 27.74%
length(cyp2c9_60x_pypgx_select)/total_cyp2c9_important_diplo

# Overal Aldy performance on important CYP2C9 alleles 
cyp2c9_10x_aldy_select <- get_diplo(correct_cyp2c9_10x_aldy, cyp2c9_important_alleles)
cyp2c9_30x_aldy_select <- get_diplo(correct_cyp2c9_30x_aldy, cyp2c9_important_alleles)
cyp2c9_60x_aldy_select <- get_diplo(correct_cyp2c9_60x_aldy, cyp2c9_important_alleles)
length(c(cyp2c9_10x_aldy_select,cyp2c9_30x_aldy_select,cyp2c9_60x_aldy_select))/
  (3*total_cyp2c9_important_diplo)
# 97.58%
length(cyp2c9_60x_aldy_select)/total_cyp2c9_important_diplo

# Overal StellarPGx performance on important CYP2C9 alleles 
cyp2c9_10x_stellarpgx_select <- get_diplo(correct_cyp2c9_10x_stellarpgx, cyp2c9_important_alleles)
cyp2c9_30x_stellarpgx_select <- get_diplo(correct_cyp2c9_30x_stellarpgx, cyp2c9_important_alleles)
cyp2c9_60x_stellarpgx_select <- get_diplo(correct_cyp2c9_60x_stellarpgx, cyp2c9_important_alleles)
length(c(cyp2c9_10x_stellarpgx_select,cyp2c9_30x_stellarpgx_select,cyp2c9_60x_stellarpgx_select))/
  (3*total_cyp2c9_important_diplo)
# 80.13%
length(cyp2c9_60x_stellarpgx_select)/total_cyp2c9_important_diplo

# For CYP2C19: *2, *3, *17
correct_cyp2c19_pypgx <- readLines("PyPGx_CYP2C19_30x_100%.txt")
correct_cyp2c19_aldy <- readLines("Aldy_CYP2C19_30x_100%.txt")
correct_cyp2c19_stellarpgx <- readLines("StellarPGx_CYP2C19_30x_100%.txt")

cyp2c19_important_alleles <- c(2, 3, 17)
total_cyp2c19_important_diplo <- sum((36-length(cyp2c19_important_alleles)+1):36)

# PyPGx performance
cyp2c19_pypgx_select <- get_diplo(correct_cyp2c19_pypgx, cyp2c19_important_alleles)
length(cyp2c19_pypgx_select)/total_cyp2c19_important_diplo
# 72.4%

# Aldy performance
cyp2c19_aldy_select <- get_diplo(correct_cyp2c19_aldy, cyp2c19_important_alleles)
length(cyp2c19_aldy_select)/total_cyp2c19_important_diplo
# 93.3%

# StellarPGx performance
cyp2c19_stellarpgx_select <- get_diplo(correct_cyp2c19_stellarpgx, cyp2c19_important_alleles)
length(cyp2c19_stellarpgx_select)/total_cyp2c19_important_diplo
# 93.3%
