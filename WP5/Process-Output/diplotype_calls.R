### Master Thesis: synthetic read simulation for PGx ###
### Lynn Henrotte 2024-2025 ###

### Packages ###

library("tidyverse")
library("ggplot2")
library("ggpubr")

###---------------------------###
###-- Find callable alleles --###
###---------------------------###

# Set working directory
setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP2")

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

###------------###
###-- CYP2C9 --###
###------------###

###-- PyPGx -- ###

### Plotting some results from the PyPGx diplotype calling for CYP2C9 ###

## Barplot indicating how correct the PyPGx diplotype calls were, at 60x coverage ##

correct <- factor(c("Wrong", "Partially Correct", "Correct"), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
number_correct <- c(3, 2122, 431)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_pypgx_calls_60x <- data.frame(correct, number_correct, percent_correct, 
                               percent_correct_label)

bar_pypgx_correct_60x <- ggplot(data = data_pypgx_calls_60x, 
                            aes(x = correct, y = percent_correct,
                                color = correct, fill = correct,
                                label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx 60x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_60x

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 3))
number_correct <- c(1, 2, 3, 2119, 67, 364)
percent_correct <- number_correct / rep(c(71, 2485), times = 3)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_pypgx_calls_60x_2 <- data.frame(correct, number_correct, percent_correct, 
                                zygosity, percent_correct_label)

bar_pypgx_correct_zygosity_60x <- ggplot(data = data_pypgx_calls_60x_2, 
                            aes(x = zygosity, y = percent_correct,
                                color = correct, fill = correct,
                                label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx 60x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_zygosity_60x

## Barplot indicating how correct the PyPGx diplotype calls were, at 30x coverage ##

correct <- factor(c("Wrong", "Partially Correct", "Correct"), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
number_correct <- c(30, 2050, 476)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_pypgx_calls_30x <- data.frame(correct, number_correct, percent_correct, 
                               percent_correct_label)

bar_pypgx_correct_30x <- ggplot(data = data_pypgx_calls_30x, 
                            aes(x = correct, y = percent_correct,
                                color = correct, fill = correct,
                                label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx 30x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_30x

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 3))
number_correct <- c(1, 29, 3, 2047, 67, 409)
percent_correct <- number_correct / rep(c(71, 2485), times = 3)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_pypgx_calls_30x_2 <- data.frame(correct, number_correct, percent_correct, 
                                     zygosity, percent_correct_label)

bar_pypgx_correct_zygosity_30x <- ggplot(data = data_pypgx_calls_30x_2, 
                                         aes(x = zygosity, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx 30x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_zygosity_30x

## Barplot indicating how correct the PyPGx diplotype calls were, at 10x coverage ##

correct <- factor(c("Wrong", "Partially Correct", "Correct"), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
number_correct <- c(5, 2015, 536)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_pypgx_calls_10x <- data.frame(correct, number_correct, percent_correct, 
                                   percent_correct_label)

bar_pypgx_correct_10x <- ggplot(data = data_pypgx_calls_10x, 
                                aes(x = correct, y = percent_correct,
                                    color = correct, fill = correct,
                                    label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx 10x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_10x

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 3))
number_correct <- c(1, 4, 0, 2015, 70, 466)
percent_correct <- number_correct / rep(c(71, 2485), times = 3)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_pypgx_calls_10x_2 <- data.frame(correct, number_correct, percent_correct, 
                                     zygosity, percent_correct_label)

bar_pypgx_correct_zygosity_10x <- ggplot(data = data_pypgx_calls_10x_2, 
                                         aes(x = zygosity, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx 10x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_zygosity_10x

## Combine zygosity-separated barplots from 10x, 30x and 60x coverage

data_pypgx_calls_10x_2$coverage <- "10x coverage"
data_pypgx_calls_30x_2$coverage <- "30x coverage"
data_pypgx_calls_60x_2$coverage <- "60x coverage"

data_pypgx_calls_ALL_2 <- rbind(data_pypgx_calls_10x_2, 
                                data_pypgx_calls_30x_2, 
                                data_pypgx_calls_60x_2)

bar_pypgx_correct_zygosity_ALL <- ggplot(data = data_pypgx_calls_ALL_2, 
                                         aes(x = zygosity, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct PyPGx CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  facet_wrap(~coverage)

bar_pypgx_correct_zygosity_ALL

## Barplot indicating how correct the PyPGx diplotype calls were, split up by coverage ##

correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), times = 3), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
coverage <- factor(rep(c("10x", "30x", "60x"), each = 3))
number_correct <- c(5, 2015, 536, 30, 2050, 476, 3, 2122, 413)
percent_correct <- number_correct / 2556
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_pypgx_calls_coverage <- data.frame(correct, number_correct, percent_correct, 
                                       coverage, percent_correct_label)

bar_pypgx_correct_coverage <- ggplot(data = data_pypgx_calls_coverage, 
                                    aes(x = coverage, y = percent_correct,
                                        color = correct, fill = correct,
                                        label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom")

bar_pypgx_correct_coverage

bar_pypgx_correct_coverage2 <- ggplot(data = data_pypgx_calls_coverage, 
                                     aes(x = correct, y = percent_correct,
                                         color = coverage, fill = coverage,
                                         label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.07) +
  labs(x = "", y = "Frequency", color = "Coverage", fill = "Coverage") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("pink3", "skyblue2", "seaolivedrab3")) +
  scale_color_manual(values = c("pink3", "skyblue2", "seaolivedrab3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_coverage2

## Barplot indicating how often each star allele was correctly called at 60x coverage

# Read in data

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
pypgx_correct_alleles_60x <- read.table("PyPGx_correct_CYP2C9_alleles_60x_freq.txt",
                                    sep = "\t", header = TRUE)

# Add zero frequency for *25
pypgx_correct_alleles_60x <- rbind(pypgx_correct_alleles_60x, c(0, 25))
pypgx_correct_alleles_60x <- pypgx_correct_alleles_60x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", pypgx_correct_alleles_60x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
pypgx_correct_alleles_60x$Allele_labels <- allele_labels

# Make color factor
allele_percent_60x <- pypgx_correct_alleles_60x$Frequency/72
pypgx_correct_alleles_60x$Percent <- allele_percent_60x
allele_colors_60x <- ifelse(allele_percent_60x > 0.85, "Good",
                        ifelse (allele_percent_60x < 0.4, "Bad", "Intermediate"))
allele_colors_60x <- factor(allele_colors_60x, 
                        levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
allele_colors2_60x <- ifelse(allele_percent_60x >= 0.5, "50% or more", "Less than 50%")
allele_colors2_60x <- factor(allele_colors2_60x, levels = c("50% or more", "Less than 50%"))

# Make a factor for the number of variants of each allele
allele_variants_60x <- c(0, rep(1, times = 70))
allele_variants_60x[c(18,35,61,68,71)] <- 2
allele_variants_60x <- factor(allele_variants_60x, levels = c("2", "1", "0"))

# Create plot with first color factor

bar_pypgx_correct_alleles_60x_1 <- ggplot(data = pypgx_correct_alleles_60x, 
                                     aes(x = allele_labels, y = Frequency,
                                         fill = allele_colors_60x,
                                         color = allele_colors_60x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 60x coverage")

bar_pypgx_correct_alleles_60x_1_percent <- ggplot(data = pypgx_correct_alleles_60x, 
                                          aes(x = allele_labels, y = Percent,
                                              fill = allele_colors_60x,
                                              color = allele_colors_60x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 60x coverage")

bar_pypgx_correct_alleles_60x_1
bar_pypgx_correct_alleles_60x_1_percent

# Create plot with second color factor

bar_pypgx_correct_alleles_60x_2 <- ggplot(data = pypgx_correct_alleles_60x, 
                                          aes(x = allele_labels, y = Frequency,
                                         fill = allele_colors2_60x,
                                         color = allele_colors2_60x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 60x coverage")

bar_pypgx_correct_alleles_60x_2_percent <- ggplot(data = pypgx_correct_alleles_60x, 
                                          aes(x = allele_labels, y = Percent,
                                              fill = allele_colors2_60x,
                                              color = allele_colors2_60x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 60x coverage")

bar_pypgx_correct_alleles_60x_2
bar_pypgx_correct_alleles_60x_2_percent

# Create plot colored by number of variants

bar_pypgx_correct_alleles_60x_3 <- ggplot(data = pypgx_correct_alleles_60x, 
                                          aes(x = allele_labels, y = Frequency,
                                         fill = allele_variants_60x,
                                         color = allele_variants_60x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Number of variants", color = "Number of variants") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("skyblue", "pink2", "grey")) +
  scale_color_manual(values = c("skyblue", "pink2", "grey")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 60x coverage")

bar_pypgx_correct_alleles_60x_3

## Barplot indicating how often each star allele was correctly called at 30x coverage

# Read in data

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
pypgx_correct_alleles_30x <- read.table("PyPGx_correct_CYP2C9_alleles_30x_freq.txt",
                                        sep = "\t", header = TRUE)

# Add zero frequency for *25
pypgx_correct_alleles_30x <- rbind(pypgx_correct_alleles_30x, c(0, 25))
pypgx_correct_alleles_30x <- pypgx_correct_alleles_30x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", pypgx_correct_alleles_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
pypgx_correct_alleles_30x$Allele_labels <- allele_labels

# Make color factor
allele_percent_30x <- pypgx_correct_alleles_30x$Frequency/72
pypgx_correct_alleles_30x$Percent <- allele_percent_30x
allele_colors_30x <- ifelse(allele_percent_30x > 0.85, "Good",
                        ifelse (allele_percent_30x < 0.4, "Bad", "Intermediate"))
allele_colors_30x <- factor(allele_colors_30x, 
                        levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
allele_colors2_30x <- ifelse(allele_percent_30x >= 0.5, "50% or more", "Less than 50%")
allele_colors2_30x <- factor(allele_colors2_30x, levels = c("50% or more", "Less than 50%"))

# Make a factor for the number of variants of each allele
allele_variants_30x <- c(0, rep(1, times = 70))
allele_variants_30x[c(18,35,61,68,71)] <- 2
allele_variants_30x <- factor(allele_variants_30x, levels = c("2", "1", "0"))

# Create plot with first color factor

bar_pypgx_correct_alleles_30x_1 <- ggplot(data = pypgx_correct_alleles_30x, 
                                          aes(x = allele_labels, y = Frequency,
                                             fill = allele_colors_30x,
                                             color = allele_colors_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", linewidth= 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_30x_1_percent <- ggplot(data = pypgx_correct_alleles_30x, 
                                          aes(x = allele_labels, y = Percent,
                                              fill = allele_colors_30x,
                                              color = allele_colors_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_30x_1
bar_pypgx_correct_alleles_30x_1_percent

# Create plot with second color factor

bar_pypgx_correct_alleles_30x_2 <- ggplot(data = pypgx_correct_alleles_30x, 
                                          aes(x = allele_labels, y = Frequency,
                                             fill = allele_colors2_30x,
                                             color = allele_colors2_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_30x_2_percent <- ggplot(data = pypgx_correct_alleles_30x, 
                                          aes(x = allele_labels, y = Percent,
                                              fill = allele_colors2_30x,
                                              color = allele_colors2_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_30x_2
bar_pypgx_correct_alleles_30x_2_percent

# Create plot colored by number of variants

bar_pypgx_correct_alleles_30x_3 <- ggplot(data = pypgx_correct_alleles_30x, 
                                          aes(x = allele_labels, y = Frequency,
                                             fill = allele_variants_30x,
                                             color = allele_variants_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Number of variants", color = "Number of variants") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("skyblue", "pink2", "grey")) +
  scale_color_manual(values = c("skyblue", "pink2", "grey")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_30x_3

## Barplot indicating how often each star allele was correctly called at 10x coverage

# Read in data

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
pypgx_correct_alleles_10x <- read.table("PyPGx_correct_CYP2C9_alleles_10x_freq.txt",
                                        sep = "\t", header = TRUE)

# Add zero frequency for *25
pypgx_correct_alleles_10x <- rbind(pypgx_correct_alleles_10x, c(0, 25))
pypgx_correct_alleles_10x <- pypgx_correct_alleles_10x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", pypgx_correct_alleles_10x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
pypgx_correct_alleles_10x$Allele_labels <- allele_labels

# Make color factor
allele_percent_10x <- pypgx_correct_alleles_10x$Frequency/72
pypgx_correct_alleles_10x$Percent <- allele_percent_10x
allele_colors_10x <- ifelse(allele_percent_10x > 0.85, "Good",
                        ifelse (allele_percent_10x < 0.4, "Bad", "Intermediate"))
allele_colors_10x <- factor(allele_colors_10x, 
                        levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
allele_colors2_10x <- ifelse(allele_percent_10x >= 0.5, "50% or more", "Less than 50%")
allele_colors2_10x <- factor(allele_colors2_10x, levels = c("50% or more", "Less than 50%"))

# Make a factor for the number of variants of each allele
allele_variants_10x <- c(0, rep(1, times = 70))
allele_variants_10x[c(18,35,61,68,71)] <- 2
allele_variants_10x <- factor(allele_variants_10x, levels = c("2", "1", "0"))

# Create plot with first color factor

bar_pypgx_correct_alleles_10x_1 <- ggplot(data = pypgx_correct_alleles_10x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_colors_10x,
                                              color = allele_colors_10x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 10x coverage")

bar_pypgx_correct_alleles_10x_1_percent <- ggplot(data = pypgx_correct_alleles_10x, 
                                          aes(x = allele_labels, y = Percent,
                                              fill = allele_colors_10x,
                                              color = allele_colors_10x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 10x coverage")

bar_pypgx_correct_alleles_10x_1
bar_pypgx_correct_alleles_10x_1_percent

# Create plot with second color factor

bar_pypgx_correct_alleles_10x_2 <- ggplot(data = pypgx_correct_alleles_10x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_colors2_10x,
                                              color = allele_colors2_10x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 10x coverage")

bar_pypgx_correct_alleles_10x_2_percent <- ggplot(data = pypgx_correct_alleles_10x, 
                                          aes(x = allele_labels, y = Percent,
                                              fill = allele_colors2_10x,
                                              color = allele_colors2_10x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Percentage called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Percentage that each CYP2C9 star allele was correctly called by PyPGx at 10x coverage")

bar_pypgx_correct_alleles_10x_2
bar_pypgx_correct_alleles_10x_2_percent

# Create plot colored by number of variants

bar_pypgx_correct_alleles_10x_3 <- ggplot(data = pypgx_correct_alleles_10x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_variants_10x,
                                              color = allele_variants_10x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 72) +
  geom_hline(yintercept = 72, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Number of variants", color = "Number of variants") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("skyblue", "pink2", "grey")) +
  scale_color_manual(values = c("skyblue", "pink2", "grey")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by PyPGx at 10x coverage")

bar_pypgx_correct_alleles_10x_3

# Combine plots for all three coverages

ggarrange(bar_pypgx_correct_alleles_10x_1,
          bar_pypgx_correct_alleles_30x_1, 
          bar_pypgx_correct_alleles_60x_1,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_pypgx_correct_alleles_10x_1_percent,
          bar_pypgx_correct_alleles_30x_1_percent, 
          bar_pypgx_correct_alleles_60x_1_percent,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_pypgx_correct_alleles_10x_2,
          bar_pypgx_correct_alleles_30x_2, 
          bar_pypgx_correct_alleles_60x_2,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_pypgx_correct_alleles_10x_2_percent,
          bar_pypgx_correct_alleles_30x_2_percent, 
          bar_pypgx_correct_alleles_60x_2_percent,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_pypgx_correct_alleles_10x_3,
          bar_pypgx_correct_alleles_30x_3, 
          bar_pypgx_correct_alleles_60x_3,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

###-- Aldy --###

## Barplot indicating how correct the Aldy diplotype calls were at 60x coverage ##

correct <- factor(c("Wrong", "Partially Correct", "Correct"), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
number_correct <- c(1, 9, 3645)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_aldy_calls_60x <- data.frame(correct, number_correct, percent_correct, 
                               percent_correct_label)

bar_aldy_correct_60x <- ggplot(data = data_aldy_calls_60x, 
                           aes(x = correct, y = percent_correct,
                               color = correct, fill = correct,
                               label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.6, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy 60x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_aldy_correct_60x

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 3))
number_correct <- c(0, 1, 9, 0, 76, 3569)
percent_correct <- number_correct / rep(c(85, 3570), times = 3)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_aldy_calls_60x_2 <- data.frame(correct, number_correct, percent_correct, 
                                     zygosity, percent_correct_label)

bar_aldy_correct_zygosity_60x <- ggplot(data = data_aldy_calls_60x_2, 
                                         aes(x = zygosity, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy 60x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_aldy_correct_zygosity_60x

## Barplot indicating how correct the Aldy diplotype calls were at 30x coverage ##

correct <- factor(c("Wrong", "Partially Correct", "Correct"), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
number_correct <- c(2, 33, 3620)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_aldy_calls_30x <- data.frame(correct, number_correct, percent_correct, 
                                  percent_correct_label)

bar_aldy_correct_30x <- ggplot(data = data_aldy_calls_30x, 
                               aes(x = correct, y = percent_correct,
                                   color = correct, fill = correct,
                                   label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.6, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy 30x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_aldy_correct_30x

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 3))
number_correct <- c(0, 2, 15, 18, 70, 3550)
percent_correct <- number_correct / rep(c(85, 3570), times = 3)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_aldy_calls_30x_2 <- data.frame(correct, number_correct, percent_correct, 
                                    zygosity, percent_correct_label)

bar_aldy_correct_zygosity_30x <- ggplot(data = data_aldy_calls_30x_2, 
                                        aes(x = zygosity, y = percent_correct,
                                            color = correct, fill = correct,
                                            label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy 30x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_aldy_correct_zygosity_30x

## Barplot indicating how correct the Aldy diplotype calls were at 10x coverage ##

correct <- factor(c("Wrong", "Partially Correct", "Correct"), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
number_correct <- c(7, 309, 3339)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_aldy_calls_10x <- data.frame(correct, number_correct, percent_correct, 
                                  percent_correct_label)

bar_aldy_correct_10x <- ggplot(data = data_aldy_calls_10x, 
                               aes(x = correct, y = percent_correct,
                                   color = correct, fill = correct,
                                   label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.6, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy 30x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_aldy_correct_10x

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 3))
number_correct <- c(1, 6, 7, 302, 77, 3262)
percent_correct <- number_correct / rep(c(85, 3570), times = 3)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_aldy_calls_10x_2 <- data.frame(correct, number_correct, percent_correct, 
                                    zygosity, percent_correct_label)

bar_aldy_correct_zygosity_10x <- ggplot(data = data_aldy_calls_10x_2, 
                                        aes(x = zygosity, y = percent_correct,
                                            color = correct, fill = correct,
                                            label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy 30x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_aldy_correct_zygosity_10x

## Barplot indicating how often each star allele was correctly called at 10x coverage

# Read in data

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
Aldy_correct_alleles_10x <- read.table("Aldy_correct_CYP2C9_alleles_10x_freq.txt",
                                             sep = "\t", header = TRUE)

# Make labels
allele_labels <- paste0("*", Aldy_correct_alleles_10x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
Aldy_correct_alleles_10x$Allele_labels <- allele_labels

# Make color factor
allele_percent_10x_Aldy <- Aldy_correct_alleles_10x$Frequency/72
allele_colors_10x_Aldy <- ifelse(allele_percent_10x_Aldy > 0.85, "Good",
                                       ifelse (allele_percent_10x_stellarpgx < 0.4, "Bad", "Intermediate"))
allele_colors_10x_stellarpgx <- factor(allele_colors_10x_stellarpgx, 
                                       levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
allele_colors2_10x_stellarpgx <- ifelse(allele_percent_10x_stellarpgx >= 0.5, "50% or more", "Less than 50%")
allele_colors2_10x_stellarpgx <- factor(allele_colors2_10x_stellarpgx, levels = c("50% or more", "Less than 50%"))

# Make a factor for the number of variants of each allele
allele_variants_10x_stellarpgx <- c(0, rep(1, times = 84))
allele_variants_10x_stellarpgx[c(18,35,61,68,71)] <- 2
allele_variants_10x_stellarpgx <- factor(allele_variants_10x_stellarpgx, levels = c("2", "1", "0"))

# Create plot with first color factor

bar_StellarPGx_correct_alleles_10x_1 <- ggplot(data = StellarPGx_correct_alleles_10x, 
                                               aes(x = allele_labels, y = Frequency,
                                                   fill = allele_colors_10x_stellarpgx,
                                                   color = allele_colors_10x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_10x_1

# Create plot with second color factor

bar_StellarPGx_correct_alleles_10x_2 <- ggplot(data = StellarPGx_correct_alleles_10x, 
                                               aes(x = allele_labels, y = Frequency,
                                                   fill = allele_colors2_10x_stellarpgx,
                                                   color = allele_colors2_10x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_10x_2

# Create plot colored by number of variants

bar_StellarPGx_correct_alleles_10x_3 <- ggplot(data = StellarPGx_correct_alleles_10x, 
                                               aes(x = allele_labels, y = Frequency,
                                                   fill = allele_variants_10x_stellarpgx,
                                                   color = allele_variants_10x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Number of variants", color = "Number of variants") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("skyblue", "pink2", "grey")) +
  scale_color_manual(values = c("skyblue", "pink2", "grey")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage")

bar_StellarPGx_correct_alleles_10x_3

## Barplot indicating how correct the Aldy diplotype calls were, split up by coverage ##

correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), times = 3), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
coverage <- factor(rep(c("10x", "30x", "60x"), each = 3))
number_correct <- c(7, 309, 3339, 2, 33, 3620, 1, 9, 3645)
percent_correct <- number_correct / 3655
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_aldy_calls_coverage <- data.frame(correct, number_correct, percent_correct, 
                                     coverage, percent_correct_label)

bar_aldy_correct_coverage <- ggplot(data = data_aldy_calls_coverage, 
                                         aes(x = coverage, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom")

bar_aldy_correct_coverage

bar_aldy_correct_coverage2 <- ggplot(data = data_aldy_calls_coverage, 
                                    aes(x = correct, y = percent_correct,
                                        color = coverage, fill = coverage,
                                        label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.07) +
  labs(x = "", y = "Frequency", color = "Coverage", fill = "Coverage") +
  ggtitle("Frequencies of wrong, partial and fully correct\nAldy CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("pink3", "skyblue2", "seaolivedrab3")) +
  scale_color_manual(values = c("pink3", "skyblue2", "seaolivedrab3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_aldy_correct_coverage2

###-- StellarPGx --###

### Plotting some results from the StellarPGx diplotype calling ###

## Barplot indicating how correct the StellarPGx diplotype calls were, at 60x coverage ##

correct <- factor(c("Failed", "Wrong", "Partially Correct", "Correct"), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
number_correct <- c(3, 75, 120, 3457)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_StellarPGx_calls_60x <- data.frame(correct, number_correct, percent_correct, 
                                   percent_correct_label)

bar_StellarPGx_correct_60x <- ggplot(data = data_StellarPGx_calls_60x, 
                                aes(x = correct, y = percent_correct,
                                    color = correct, fill = correct,
                                    label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx 60x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_StellarPGx_correct_60x

## Ideas for further enhancement of this plot:
## - split up bars by SNPs/indels
## - split up bars by simulation method (if others are used)

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Failed", "Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 4))
number_correct <- c(0, 3, 1, 74, 0, 120, 84, 3373)
percent_correct <- number_correct / rep(c(85, 3570), times = 4)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_StellarPGx_calls_60x_2 <- data.frame(correct, number_correct, percent_correct, 
                                     zygosity, percent_correct_label)

bar_StellarPGx_correct_zygosity_60x <- ggplot(data = data_StellarPGx_calls_60x_2, 
                                         aes(x = zygosity, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx 60x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_StellarPGx_correct_zygosity_60x

## Barplot indicating how correct the StellarPGx diplotype calls were, at 30x coverage ##

correct <- factor(c("Failed", "Wrong", "Partially Correct", "Correct"), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
number_correct <- c(6, 38, 617, 2994)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_StellarPGx_calls_30x <- data.frame(correct, number_correct, percent_correct, 
                                   percent_correct_label)

bar_StellarPGx_correct_30x <- ggplot(data = data_StellarPGx_calls_30x, 
                                aes(x = correct, y = percent_correct,
                                    color = correct, fill = correct,
                                    label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx 30x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_StellarPGx_correct_30x

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Failed", "Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 4))
number_correct <- c(0, 6, 1, 37, 0, 617, 84, 2910)
percent_correct <- number_correct / rep(c(85, 3570), times = 4)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_StellarPGx_calls_30x_2 <- data.frame(correct, number_correct, percent_correct, 
                                     zygosity, percent_correct_label)

bar_StellarPGx_correct_zygosity_30x <- ggplot(data = data_StellarPGx_calls_30x_2, 
                                         aes(x = zygosity, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx 30x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_StellarPGx_correct_zygosity_30x

## Barplot indicating how correct the StellarPGx diplotype calls were, at 10x coverage ##

correct <- factor(c("Failed", "Wrong", "Partially Correct", "Correct"), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
number_correct <- c(5, 148, 1288, 2214)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

data_StellarPGx_calls_10x <- data.frame(correct, number_correct, percent_correct, 
                                   percent_correct_label)

bar_StellarPGx_correct_10x <- ggplot(data = data_StellarPGx_calls_10x, 
                                aes(x = correct, y = percent_correct,
                                    color = correct, fill = correct,
                                    label = percent_correct_label)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) + theme_bw() + ylim(0,1.05) +
  guides(color = "none", fill = "none") +
  labs(x = "", y = "Frequency") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx 10x CYP2C9 diplotype calls") +
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5))

bar_StellarPGx_correct_10x

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Failed", "Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 4))
number_correct <- c(0, 5, 2, 146, 2, 1286, 81, 2133)
percent_correct <- number_correct / rep(c(85, 3570), times = 4)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_StellarPGx_calls_10x_2 <- data.frame(correct, number_correct, percent_correct, 
                                     zygosity, percent_correct_label)

bar_StellarPGx_correct_zygosity_10x <- ggplot(data = data_StellarPGx_calls_10x_2, 
                                         aes(x = zygosity, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx 10x CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_StellarPGx_correct_zygosity_10x

## Combine zygosity-separated barplots from 10x, 30x and 60x coverage

data_StellarPGx_calls_10x_2$coverage <- "10x coverage"
data_StellarPGx_calls_30x_2$coverage <- "30x coverage"
data_StellarPGx_calls_60x_2$coverage <- "60x coverage"

data_StellarPGx_calls_ALL_2 <- rbind(data_StellarPGx_calls_10x_2, 
                                data_StellarPGx_calls_30x_2, 
                                data_StellarPGx_calls_60x_2)

bar_StellarPGx_correct_zygosity_ALL <- ggplot(data = data_StellarPGx_calls_ALL_2, 
                                         aes(x = zygosity, y = percent_correct,
                                             color = correct, fill = correct,
                                             label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx CYP2C9 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  facet_wrap(~coverage)

bar_StellarPGx_correct_zygosity_ALL

## Barplot indicating how correct the StellarPGx diplotype calls were, split up by coverage ##

correct <- factor(rep(c("Failed", "Wrong", "Partially Correct", "Correct"), times = 3), 
                  levels = c("Failed", "Wrong", "Partially Correct", "Correct"))
coverage <- factor(rep(c("10x", "30x", "60x"), each = 4))
number_correct <- c(5, 148, 1288, 2214, 6, 38, 617, 2994, 3, 75, 120, 3457)
percent_correct <- number_correct / 3655
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_StellarPGx_calls_coverage <- data.frame(correct, number_correct, percent_correct, 
                                        coverage, percent_correct_label)

bar_StellarPGx_correct_coverage <- ggplot(data = data_StellarPGx_calls_coverage, 
                                     aes(x = coverage, y = percent_correct,
                                         color = correct, fill = correct,
                                         label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("black", "firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        legend.position = "bottom")

bar_StellarPGx_correct_coverage

bar_StellarPGx_correct_coverage2 <- ggplot(data = data_StellarPGx_calls_coverage, 
                                      aes(x = correct, y = percent_correct,
                                          color = coverage, fill = coverage,
                                          label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.07) +
  labs(x = "", y = "Frequency", color = "Coverage", fill = "Coverage") +
  ggtitle("Frequencies of wrong, partial and fully correct\nStellarPGx CYP2C9 diplotype calls, by coverage") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("pink3", "skyblue2", "seagreen3")) +
  scale_color_manual(values = c("pink3", "skyblue2", "seagreen3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_StellarPGx_correct_coverage2

### CONCLUSIONS BASED ON THIS GRAPH: AS COVERAGE DECREASES, THE NUMBER OF CORRECT DIPLOTYPE CALLS BY STELLARPGX
### DECREASES AS WELL. HOWEVER, IT IS PROMISING THAT MOST OF THESE MISSED CALLS END UP AS PARTIALLY CORRECT
### INSTEAD OF COMPLETELY WRONG, THE CATEGORY WHICH DOES NOT SEEM TO CHANGE A LOT WHEN COVERAGE CHANGES.


## Barplot indicating how often each star allele was correctly called at 60x coverage

# Read in data

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
StellarPGx_correct_alleles_60x <- read.table("StellarPGx_correct_CYP2C9_alleles_60x_freq.txt",
                                        sep = "\t", header = TRUE)

# Add zero frequency for *25
StellarPGx_correct_alleles_60x <- rbind(StellarPGx_correct_alleles_60x, c(0, 25))
StellarPGx_correct_alleles_60x <- StellarPGx_correct_alleles_60x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", StellarPGx_correct_alleles_60x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
StellarPGx_correct_alleles_60x$Allele_labels <- allele_labels

# Make color factor
allele_percent_60x_stellarpgx <- StellarPGx_correct_alleles_60x$Frequency/86
allele_colors_60x_stellarpgx <- ifelse(allele_percent_60x_stellarpgx > 0.85, "Good",
                            ifelse (allele_percent_60x_stellarpgx < 0.4, "Bad", "Intermediate"))
allele_colors_60x_stellarpgx <- factor(allele_colors_60x_stellarpgx, 
                            levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
allele_colors2_60x_stellarpgx <- ifelse(allele_percent_60x_stellarpgx >= 0.5, "50% or more", "Less than 50%")
allele_colors2_60x_stellarpgx <- factor(allele_colors2_60x_stellarpgx, levels = c("50% or more", "Less than 50%"))

# Make a factor for the number of variants of each allele
allele_variants_60x_stellarpgx <- c(0, rep(1, times = 84))
allele_variants_60x_stellarpgx[c(18,35,61,68,71)] <- 2
allele_variants_60x_stellarpgx <- factor(allele_variants_60x_stellarpgx, levels = c("2", "1", "0"))

# Create plot with first color factor

bar_StellarPGx_correct_alleles_60x_1 <- ggplot(data = StellarPGx_correct_alleles_60x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_colors_60x_stellarpgx,
                                              color = allele_colors_60x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 60x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_60x_1

# Create plot with second color factor

bar_StellarPGx_correct_alleles_60x_2 <- ggplot(data = StellarPGx_correct_alleles_60x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_colors2_60x_stellarpgx,
                                              color = allele_colors2_60x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 60x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_60x_2

# Create plot colored by number of variants

bar_StellarPGx_correct_alleles_60x_3 <- ggplot(data = StellarPGx_correct_alleles_60x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_variants_60x_stellarpgx,
                                              color = allele_variants_60x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Number of variants", color = "Number of variants") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("skyblue", "pink2", "grey")) +
  scale_color_manual(values = c("skyblue", "pink2", "grey")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 60x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_60x_3

## Barplot indicating how often each star allele was correctly called at 30x coverage

# Read in data

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
StellarPGx_correct_alleles_30x <- read.table("StellarPGx_correct_CYP2C9_alleles_30x_freq.txt",
                                        sep = "\t", header = TRUE)

# Make labels
allele_labels <- paste0("*", StellarPGx_correct_alleles_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
StellarPGx_correct_alleles_30x$Allele_labels <- allele_labels

# Make color factor
allele_percent_30x_stellarpgx <- StellarPGx_correct_alleles_30x$Frequency/86
allele_colors_30x_stellarpgx <- ifelse(allele_percent_30x_stellarpgx > 0.85, "Good",
                            ifelse (allele_percent_30x_stellarpgx < 0.4, "Bad", "Intermediate"))
allele_colors_30x_stellarpgx <- factor(allele_colors_30x_stellarpgx, 
                            levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
allele_colors2_30x_stellarpgx <- ifelse(allele_percent_30x_stellarpgx >= 0.5, "50% or more", "Less than 50%")
allele_colors2_30x_stellarpgx <- factor(allele_colors2_30x_stellarpgx, levels = c("50% or more", "Less than 50%"))

# Make a factor for the number of variants of each allele
allele_variants_30x_stellarpgx <- c(0, rep(1, times = 84))
allele_variants_30x_stellarpgx[c(18,35,61,68,71)] <- 2
allele_variants_30x_stellarpgx <- factor(allele_variants_30x_stellarpgx, levels = c("2", "1", "0"))

# Create plot with first color factor

bar_StellarPGx_correct_alleles_30x_1 <- ggplot(data = StellarPGx_correct_alleles_30x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_colors_30x_stellarpgx,
                                              color = allele_colors_30x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 30x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_30x_1

# Create plot with second color factor

bar_StellarPGx_correct_alleles_30x_2 <- ggplot(data = StellarPGx_correct_alleles_30x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_colors2_30x_stellarpgx,
                                              color = allele_colors2_30x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 30x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_30x_2

# Create plot colored by number of variants

bar_StellarPGx_correct_alleles_30x_3 <- ggplot(data = StellarPGx_correct_alleles_30x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_variants_30x_stellarpgx,
                                              color = allele_variants_30x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Number of variants", color = "Number of variants") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("skyblue", "pink2", "grey")) +
  scale_color_manual(values = c("skyblue", "pink2", "grey")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 30x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_30x_3

## Barplot indicating how often each star allele was correctly called at 10x coverage

# Read in data

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
StellarPGx_correct_alleles_10x <- read.table("StellarPGx_correct_CYP2C9_alleles_10x_freq.txt",
                                        sep = "\t", header = TRUE)

# Add zero frequency for *25
StellarPGx_correct_alleles_10x <- rbind(StellarPGx_correct_alleles_10x, c(0, 25))
StellarPGx_correct_alleles_10x <- StellarPGx_correct_alleles_10x %>% arrange(Allele)

# Make labels
allele_labels <- paste0("*", StellarPGx_correct_alleles_10x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
StellarPGx_correct_alleles_10x$Allele_labels <- allele_labels

# Make color factor
allele_percent_10x_stellarpgx <- StellarPGx_correct_alleles_10x$Frequency/72
allele_colors_10x_stellarpgx <- ifelse(allele_percent_10x_stellarpgx > 0.85, "Good",
                            ifelse (allele_percent_10x_stellarpgx < 0.4, "Bad", "Intermediate"))
allele_colors_10x_stellarpgx <- factor(allele_colors_10x_stellarpgx, 
                            levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
allele_colors2_10x_stellarpgx <- ifelse(allele_percent_10x_stellarpgx >= 0.5, "50% or more", "Less than 50%")
allele_colors2_10x_stellarpgx <- factor(allele_colors2_10x_stellarpgx, levels = c("50% or more", "Less than 50%"))

# Make a factor for the number of variants of each allele
allele_variants_10x_stellarpgx <- c(0, rep(1, times = 84))
allele_variants_10x_stellarpgx[c(18,35,61,68,71)] <- 2
allele_variants_10x_stellarpgx <- factor(allele_variants_10x_stellarpgx, levels = c("2", "1", "0"))

# Create plot with first color factor

bar_StellarPGx_correct_alleles_10x_1 <- ggplot(data = StellarPGx_correct_alleles_10x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_colors_10x_stellarpgx,
                                              color = allele_colors_10x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_10x_1

# Create plot with second color factor

bar_StellarPGx_correct_alleles_10x_2 <- ggplot(data = StellarPGx_correct_alleles_10x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_colors2_10x_stellarpgx,
                                              color = allele_colors2_10x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage")

#### ADD PERCENTAGES ON Y AXIS

bar_StellarPGx_correct_alleles_10x_2

# Create plot colored by number of variants

bar_StellarPGx_correct_alleles_10x_3 <- ggplot(data = StellarPGx_correct_alleles_10x, 
                                          aes(x = allele_labels, y = Frequency,
                                              fill = allele_variants_10x_stellarpgx,
                                              color = allele_variants_10x_stellarpgx)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 86) +
  geom_hline(yintercept = 86, color = "tomato", size = 1, lty = 2) +
  labs(x = "CYP2C9 star allele", y = "Times called", 
       fill = "Number of variants", color = "Number of variants") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("skyblue", "pink2", "grey")) +
  scale_color_manual(values = c("skyblue", "pink2", "grey")) +
  ggtitle("Number of times each CYP2C9 star allele was correctly called by StellarPGx at 10x coverage")

bar_StellarPGx_correct_alleles_10x_3

# Combine plots for all three coverages

ggarrange(bar_StellarPGx_correct_alleles_10x_1,
          bar_StellarPGx_correct_alleles_30x_1, 
          bar_StellarPGx_correct_alleles_60x_1,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_StellarPGx_correct_alleles_10x_2,
          bar_StellarPGx_correct_alleles_30x_2, 
          bar_StellarPGx_correct_alleles_60x_2,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

ggarrange(bar_StellarPGx_correct_alleles_10x_3,
          bar_StellarPGx_correct_alleles_30x_3, 
          bar_StellarPGx_correct_alleles_60x_3,
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")

## Conclusions based on these plots:
## As coverage decreases, more star alleles are more commonly not correctly called.
## However, despite *25, there does not seem to be any more alleles that cause
## StellarPGx trouble in general (i.e. only for one of the coverages)

###-- Comparisons --###

### First, compare barplots with correctness of diplotype calls

ggarrange(bar_pypgx_correct_60x, bar_aldy_correct_60x, ncol = 2, nrow = 1)
ggarrange(bar_pypgx_correct_30x, bar_aldy_correct_30x, ncol = 2, nrow = 1)

ggarrange(bar_pypgx_correct_60x, bar_aldy_correct_60x, 
          bar_pypgx_correct_30x, bar_aldy_correct_30x, ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom")

### Next, compare barplots with correctness of diplotype calls, split up by coverage

ggarrange(bar_pypgx_correct_coverage, bar_aldy_correct_coverage, nrow = 1,
          common.legend = TRUE, legend = "bottom")

ggarrange(bar_pypgx_correct_coverage2, bar_aldy_correct_coverage2, nrow = 1,
          common.legend = TRUE, legend = "bottom")

### Compare barplots with correctness of diplotype calls, split up by zygosity

ggarrange(bar_pypgx_correct_zygosity_60x, bar_aldy_correct_zygosity_60x,
          bar_StellarPGx_correct_zygosity_60x, 
          ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")

ggarrange(bar_pypgx_correct_zygosity_30x, 
          bar_aldy_correct_zygosity_30x,
          bar_StellarPGx_correct_zygosity_30x,
          common.legend = TRUE, legend = "bottom")

ggarrange(bar_pypgx_correct_zygosity_60x, bar_aldy_correct_zygosity_60x, bar_StellarPGx_correct_zygosity_60x,
          bar_pypgx_correct_zygosity_30x, bar_aldy_correct_zygosity_30x, bar_StellarPGx_correct_zygosity_30x,
          bar_pypgx_correct_zygosity_10x, bar_aldy_correct_zygosity_10x, bar_StellarPGx_correct_zygosity_10x,
          common.legend = TRUE, legend = "bottom")

#------------------------------------------------------------------------------------------------------#

###-------------###
###-- CYP2C19 --###
###-------------###

###-- PyPGx --###

## Barplot indicating how correct the PyPGx diplotype calls were, at 30x coverage ##

correct <- factor(c("Wrong", "Partially Correct", "Correct"), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
number_correct <- c(3, 296, 367)
percent_correct <- number_correct / sum(number_correct)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "% (", number_correct, ")")

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
  geom_text(position = position_dodge(width = 0.8), vjust = -0.5, size = 8,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_30x_cyp2c19

## Same plot, but split up bars by homo/heterozygosity of the diplotypes
correct <- factor(rep(c("Wrong", "Partially Correct", "Correct"), each = 2), 
                  levels = c("Wrong", "Partially Correct", "Correct"))
zygosity <- factor(rep(c("Homozygous", "Heterozygous"), times = 3))
number_correct <- c(2, 1, 0, 296, 34, 333)
percent_correct <- number_correct / rep(c(36, 630), times = 3)
percent_correct_label <- paste0(round(percent_correct * 100, 2), 
                                "%\n(", number_correct, ")")

data_pypgx_calls_30x_cyp2c19_2 <- data.frame(correct, number_correct, percent_correct, 
                                             zygosity, percent_correct_label)

bar_pypgx_correct_zygosity_30x_cyp2c19 <- ggplot(data = data_pypgx_calls_30x_cyp2c19_2, 
                                                 aes(x = zygosity, y = percent_correct,
                                                     color = correct, fill = correct,
                                                     label = percent_correct_label)) +
  geom_bar(position = "dodge", stat = "identity", alpha = 0.7, width = 0.8) + 
  theme_bw() + ylim(0,1.05) +
  labs(x = "", y = "Frequency", color = "Call status", fill = "Call status") +
  ggtitle("Frequencies of wrong, partial and fully correct\nPyPGx 30x CYP2C19 diplotype calls, by zygosity") +
  geom_text(position = position_dodge(width = 0.8), vjust = -.2, size = 5,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = c("firebrick3", "orange2", "green3")) +
  scale_color_manual(values = c("firebrick3", "orange2", "green3")) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

bar_pypgx_correct_zygosity_30x_cyp2c19

## Barplot indicating how often each star allele was correctly called at 30x coverage

# Read in data

setwd("C:/Users/lynnh/OneDrive/Bureaublad/2nd Master Stat/Master Thesis/WP5")
pypgx_correct_alleles_cyp2c19_30x <- read.table("PyPGx_correct_CYP2C19_alleles_30x_freq.txt",
                                                sep = "\t", header = TRUE)
pypgx_correct_alleles_cyp2c19_30x <- arrange(pypgx_correct_alleles_cyp2c19_30x, Allele)

# Adjust frequencies 
new_freq <- pypgx_correct_alleles_cyp2c19_30x$Frequency
new_freq <- ifelse(new_freq > 0, new_freq - 5, new_freq)
pypgx_correct_alleles_cyp2c19_30x$Frequency <- new_freq

# Make labels
allele_labels <- paste0("*", pypgx_correct_alleles_cyp2c19_30x$Allele)
allele_labels <- factor(allele_labels, levels = allele_labels)
pypgx_correct_alleles_cyp2c19_30x$Allele_labels <- allele_labels

# Make color factor
allele_percent_cyp2c19_30x <- pypgx_correct_alleles_cyp2c19_30x$Frequency/37
pypgx_correct_alleles_cyp2c19_30x$Percent <- allele_percent_cyp2c19_30x
allele_colors_cyp2c19_30x <- ifelse(allele_percent_cyp2c19_30x > 0.85, "Good",
                                    ifelse (allele_percent_cyp2c19_30x < 0.4, "Bad", "Intermediate"))
allele_colors_cyp2c19_30x <- factor(allele_colors_cyp2c19_30x, 
                                    levels = c("Good", "Intermediate", "Bad"))

# Make another color factor
allele_colors2_cyp2c19_30x <- ifelse(allele_percent_cyp2c19_30x >= 0.5, "50% or more", "Less than 50%")
allele_colors2_cyp2c19_30x <- factor(allele_colors2_cyp2c19_30x, levels = c("50% or more", "Less than 50%"))

# Make a factor for the number of variants of each allele
allele_variants_cyp2c19_30x <- c(rep(c(2), times = 32), "SV", "SV", 0, 3)
allele_variants_cyp2c19_30x[c(1,16,27)] <- 1
allele_variants_cyp2c19_30x[c(2,25)] <- 3
allele_variants_cyp2c19_30x <- factor(allele_variants_cyp2c19_30x, levels = c("SV", "3", "2", "1", "0"))

# Create plot with first color factor

bar_pypgx_correct_alleles_cyp2c19_30x_1 <- ggplot(data = pypgx_correct_alleles_cyp2c19_30x, 
                                                  aes(x = allele_labels, y = Frequency,
                                                      fill = allele_colors_cyp2c19_30x,
                                                      color = allele_colors_cyp2c19_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 37) +
  geom_hline(yintercept = 37, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Times called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C19 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_cyp2c19_30x_1_percent <- ggplot(data = pypgx_correct_alleles_cyp2c19_30x, 
                                                          aes(x = allele_labels, y = Percent,
                                                              fill = allele_colors_cyp2c19_30x,
                                                              color = allele_colors_cyp2c19_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C19 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_cyp2c19_30x_1
bar_pypgx_correct_alleles_cyp2c19_30x_1_percent

# Create plot with second color factor

bar_pypgx_correct_alleles_cyp2c19_30x_2 <- ggplot(data = pypgx_correct_alleles_cyp2c19_30x, 
                                                  aes(x = allele_labels, y = Frequency,
                                                      fill = allele_colors2_cyp2c19_30x,
                                                      color = allele_colors2_cyp2c19_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 37) +
  geom_hline(yintercept = 37, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Times called", 
       fill = "Percent allele called", color = "Percent allele called") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab3", "firebrick3")) +
  scale_color_manual(values = c("olivedrab3", "firebrick3")) +
  ggtitle("Number of times each CYP2C19 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_cyp2c19_30x_2_percent <- ggplot(data = pypgx_correct_alleles_cyp2c19_30x, 
                                                          aes(x = allele_labels, y = Percent,
                                                              fill = allele_colors_cyp2c19_30x,
                                                              color = allele_colors_cyp2c19_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 1) +
  geom_hline(yintercept = 1, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Percentage called", 
       fill = "Allele call status", color = "Allele call status") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                    values = c("olivedrab3", "grey", "firebrick3")) +
  scale_color_manual(labels = c("Good (>85%)", "Intermediate", "Bad (<40%)"),
                     values = c("olivedrab3", "grey", "firebrick3")) +
  ggtitle("Number of times each CYP2C19 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_cyp2c19_30x_2
bar_pypgx_correct_alleles_cyp2c19_30x_2_percent

# Create plot colored by number of variants

bar_pypgx_correct_alleles_cyp2c19_30x_3 <- ggplot(data = pypgx_correct_alleles_cyp2c19_30x, 
                                                  aes(x = allele_labels, y = Frequency,
                                                      fill = allele_variants_cyp2c19_30x,
                                                      color = allele_variants_cyp2c19_30x)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.7) + 
  theme_bw() + ylim(0, 37) +
  geom_hline(yintercept = 37, color = "tomato", linewidth = 1, lty = 2) +
  labs(x = "CYP2C19 star allele", y = "Times called", 
       fill = "Number of variants", color = "Number of variants") + 
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5),
        text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  scale_fill_manual(values = c("olivedrab", "skyblue", "pink2", "tomato3", "grey")) +
  scale_color_manual(values = c("olivedrab", "skyblue", "pink2", "tomato3", "grey")) +
  ggtitle("Number of times each CYP2C19 star allele was correctly called by PyPGx at 30x coverage")

bar_pypgx_correct_alleles_cyp2c19_30x_3
