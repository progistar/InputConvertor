install.packages("ggVennDiagram")
library(ggVennDiagram)
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

get_peptide_seq <- function(data, type) {
  data <- data[data["IsCanonical"] == type, ]
  data <- data[!duplicated(data$InferredPeptide), ]
  return (data$InferredPeptide)
}

get_peptide_len <- function(data, type, sample) {
  data <- data[data["IsCanonical"] == type, ]
  data <- data[!duplicated(data$InferredPeptide), ]
  data$len <- nchar(data$InferredPeptide)
  data <- data.frame(
    values = c(nrow(data[data$len==7,]), 
               nrow(data[data$len==8,]),
               nrow(data[data$len==9,]),
               nrow(data[data$len==10,]),
               nrow(data[data$len==11,]),
               nrow(data[data$len==12,]),
               nrow(data[data$len==13,]),
               nrow(data[data$len==14,]),
               nrow(data[data$len==15,]))
  )
  
  data <- t(data)
  colnames(data) <- c(7,8,9,10,11,12,13,14,15)
  
  return (data)
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
  
  
  ######### Peptide overlap ##############
  # Create a list with the elements of each set
  set1 <- get_peptide_seq(CSNV, "true")
  set2 <- get_peptide_seq(PNOVO3, "true")
  set3 <- get_peptide_seq(PEAKS, "true")
  # Create the Venn diagram
  x <- list(Casanovo = set1, pNovo3 = set2, PEAKS11 = set3)
  plot <- ggVennDiagram(x) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
    ggtitle(paste("Canonical peptides in ", sample, sep = "")) + theme_minimal() + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  ggsave(paste("cPeptsIn", sample, ".png", sep = ""), plot, width = 8, height = 6, dpi = 600)
  
  set1 <- get_peptide_seq(CSNV, "false")
  set2 <- get_peptide_seq(PNOVO3, "false")
  set3 <- get_peptide_seq(PEAKS, "false")
  # Create the Venn diagram
  x <- list(Casanovo = set1, pNovo3 = set2, PEAKS11 = set3)
  plot <- ggVennDiagram(x) + scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") + 
    ggtitle(paste("Noncanonical peptides in ", sample, sep = "")) + theme_minimal() + theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )
  ggsave(paste("ncPeptsIn", sample, ".png", sep = ""), plot, width = 8, height = 6, dpi = 600)
  
  
  ########## Peptide length distribution ###########
  set1 <- get_peptide_len(CSNV, "true")
  set2 <- get_peptide_len(PNOVO3, "true")
  set3 <- get_peptide_len(PEAKS, "true")
  
  len <- c(7,8,9,10,11,12,13,14,15)
  df <- data.frame(
    Value = c(set1, set2, set3),
    Length = c(len, len, len),
    Sample = rep(c("Casanovo", "pNovo3", "PEAKS11"), each = 9)
  )
  
  df$Sample <- factor(df$Sample, levels = c("Casanovo", "pNovo3", "PEAKS11"))
  
  plot <- ggplot(df, aes(x = factor(Length), y= Value, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    labs(title = paste("Canonical peptides in ", sample, sep = ""), x = "Length", y = "No. peptides") +
    scale_x_discrete(label = len) +
    scale_fill_manual(values = c("Casanovo" = "skyblue", "pNovo3" = "lightgreen", "PEAKS11" = "lightcoral")) + theme_minimal()
  ggsave(paste("cPeptsLenIn", sample, ".png", sep = ""), plot, width = 8, height = 6, dpi = 600)
  
  set1 <- get_peptide_len(CSNV, "false")
  set2 <- get_peptide_len(PNOVO3, "false")
  set3 <- get_peptide_len(PEAKS, "false")
  
  len <- c(7,8,9,10,11,12,13,14,15)
  df <- data.frame(
    Value = c(set1, set2, set3),
    Length = c(len, len, len),
    Sample = rep(c("Casanovo", "pNovo3", "PEAKS11"), each = 9)
  )
  
  df$Sample <- factor(df$Sample, levels = c("Casanovo", "pNovo3", "PEAKS11"))
  
  plot <- ggplot(df, aes(x = factor(Length), y= Value, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    labs(title = paste("Noncanonical peptides in ", sample, sep = ""), x = "Length", y = "No. peptides") +
    scale_x_discrete(label = len) +
    scale_fill_manual(values = c("Casanovo" = "skyblue", "pNovo3" = "lightgreen", "PEAKS11" = "lightcoral")) + theme_minimal()
  ggsave(paste("ncPeptsLenIn", sample, ".png", sep = ""), plot, width = 8, height = 6, dpi = 600)
  
}


colnames(data_summary) <- c("Sample", "Tool", "Class", "No. PSMs", "Binder ratio")
colnames(data_summary)

data_summary$`No. PSMs` <- as.numeric(data_summary$`No. PSMs`)
data_summary$`Binder ratio` <- as.numeric(data_summary$`Binder ratio`)

plot <- ggplot(data_summary[data_summary$Class == "Canonical",], 
       aes(x = `Binder ratio`, y = `No. PSMs`, shape = factor(Tool), color = factor(Sample))) +
  geom_point(size = 3) +
  scale_shape_manual(name = "Tool", values = c(1,2,3)) +
  scale_color_manual(name = "Sample", values = c("red", "blue", "green", "black")) +
  labs(title = "Performance test for canonical PSMs") +
  theme_minimal()

ggsave("performance_canonical.png", plot, width = 8, height = 6, dpi = 600)




