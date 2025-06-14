library("tidyverse")

data <- read.delim("guidelineAnnotationTable-all-data.tsv")

data_dpwg <- data %>% filter(dpwg != "")

data_common <- data %>% filter(cpic != "") %>% filter(dpwg != "")
length(unique(data_common$drug))
