library(ggplot2)

get_number_of_psms <- function(data, type) {
  data <- data[data["IsCanonical"] == type, ]
  return (nrow(data))
}

get_ba_ratio_of_psms <- function(data, type) {
  data <- data[data["IsCanonical"] == type, ]
  data <- data[data$BestScore != "-", ] ## discard out of lengths
  
  
  
  return (nrow(data[data$HLA != "NB", ]) / nrow(data))
}

setwd("/Users/gistar/Documents/ZhangLab/2023_Immunopeptidomics_LUAD/Test/4sample_test_fdr")

samples <- c("KZ2022_C3N-01488.T", "KZ2022_C3N-02423.T",
             "VY2021_C3L-02549.T", "VY2021_C3N-01016.T")

data_summary <- as.data.frame(matrix(ncol=5, nrow=0), stringsAsFactors = FALSE, )
for(sample in samples) {
  print(sample)
  
  CSNV <- read.csv(
    paste(sample, ".csnv.pXg.feat.fdr.ba", sep = ""),
    header = T,
    sep="\t"
  )
  
  PEAKS <- read.csv(
    paste(sample, ".PEAKS11.pXg.feat.fdr.ba", sep = ""),
    header = T,
    sep="\t"
  )
  
  PNOVO3 <- read.csv(
    paste(sample, ".pNovo3.pXg.feat.fdr.ba", sep = ""),
    header = T,
    sep="\t"
  )
  
  data_summary <- rbind(data_summary, 
                        c(sample,"Casanovo", "Canonical", get_number_of_psms(CSNV, "true"), get_ba_ratio_of_psms(CSNV, "true")))
  data_summary <- rbind(data_summary, 
                        c(sample,"Casanovo", "Noncanonical", get_number_of_psms(CSNV, "false"), get_ba_ratio_of_psms(CSNV, "false")))
  data_summary <- rbind(data_summary, 
                        c(sample,"PEAKS11", "Canonical", get_number_of_psms(PEAKS, "true"), get_ba_ratio_of_psms(PEAKS, "true")))
  data_summary <- rbind(data_summary, 
                        c(sample,"PEAKS11", "Noncanonical", get_number_of_psms(PEAKS, "false"), get_ba_ratio_of_psms(PEAKS, "false")))
  data_summary <- rbind(data_summary, 
                        c(sample,"pNovo3", "Canonical", get_number_of_psms(PNOVO3, "true"), get_ba_ratio_of_psms(PNOVO3, "true")))
  data_summary <- rbind(data_summary, 
                        c(sample,"pNovo3", "Noncanonical", get_number_of_psms(PNOVO3, "false"), get_ba_ratio_of_psms(PNOVO3, "false")))
  
}

colnames(data_summary) <- c("Sample", "Tool", "Class", "No. PSMs", "Binder ratio")
colnames(data_summary)

data_summary$`No. PSMs` <- as.numeric(data_summary$`No. PSMs`)
data_summary$`Binder ratio` <- as.numeric(data_summary$`Binder ratio`)

ggplot(data_summary[data_summary$Class == "Noncanonical",], 
       aes(x = `Binder ratio`, y = `No. PSMs`, shape = factor(Tool), color = factor(Sample))) +
  geom_point(size = 3) +
  scale_shape_manual(name = "Tool", values = c(1,2,3)) +
  scale_color_manual(name = "Sample", values = c("red", "blue", "green", "black")) +
  labs(title = "Performance test for noncanonical PSMs") +
  theme_minimal()





