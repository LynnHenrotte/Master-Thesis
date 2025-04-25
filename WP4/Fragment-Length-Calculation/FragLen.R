###-------------------------------###
### Master Thesis - Lynn Henrotte ###
###-------------------------------###############################
### Estimation of fragment length mean and standard deviation ###
###-----------------------------------------------------------###

## Packages ##

library("tidyverse")

## Preliminaries ##
path_to_folder = paste0("C:/Users/lynnh/OneDrive/Bureaublad",
                        "/2nd Master Stat/Master Thesis/WP4/",
                        "Fragment-Length-Calculation")
setwd(path_to_folder)

## Obtain data ##

hist_files <- list.files(path = path_to_folder, pattern = "*.csv$")

samples_list <- list()
sample_names <- c()
for (i in 1:length(hist_files)) {
  filename <- hist_files[i]
  data <- read.csv(filename, skip = 1, header = TRUE)
  colnames(data)[2] <- paste0("Count",i)
  samples_list[[i]] <- data
  sample_names[i] <- paste0("sample",i)
}
names(samples_list) <- sample_names
rm(data)

## Merge samples

samples_merged <- samples_list[[1]]

for (i in 2:length(samples_list)) {
  samples_merged <- merge(samples_merged, samples_list[[i]], 
                          id = "FragmentLength", all = TRUE)
}

samples_merged[is.na(samples_merged)] <- 0

samples_merged$Count <- rowSums(samples_merged[,-c(1)])

data_final <- samples_merged %>% 
  dplyr::select(FragmentLength, Count) %>%
  filter(Count > 0)

rm(samples_merged)

## Compute mean fragment length

data_final$Prod <- data_final$FragmentLength * data_final$Count
TotalNumFragments <- sum(data_final$Count)
MeanFragLen <- sum(data_final$Prod)/TotalNumFragments

## Compute standard deviation of fragment length

data_final$Resid <- (data_final$FragmentLength - MeanFragLen)^2
data_final$Resid_Prod <- data_final$Resid * data_final$Count
VarFragLen <- sum(data_final$Resid_Prod)/(TotalNumFragments - 1)
sqrt(VarFragLen)

## Remove fragments with length > 800
## Justification: any fragment of length greater than 800 is likely false, due
## to alignment issues.
data_crop <- data_final %>% filter(FragmentLength <= 800)
var_crop <- sum(data_crop$Resid_Prod)/(sum(data_crop$Count) - 1)
sqrt(var_crop)
