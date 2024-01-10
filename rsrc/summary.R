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

get_ba_ratio_of_peptides <- function(data, type) {
  data <- data[data["IsCanonical"] == type, ]
  data <- data[data$BestScore != "-", ] ## discard out of lengths
  
  total_pept <- get_peptide_seq(data, type)
  ba_pept <- get_peptide_seq(data[data$HLA != "NB", ], type)
  
  
  return (length(ba_pept) / length(total_pept))
}

get_peptide_seq <- function(data, type) {
  data <- data[data["IsCanonical"] == type, ]
  data <- data[!duplicated(data$InferredPeptide), ]
  return (data$InferredPeptide)
}

get_number_of_peptides <- function(data, type) {
  data <- get_peptide_seq(data, type)
  return (length(data))
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

setwd("/Users/seunghyukchoi/Documents/1_Projects/2023_Neoflow2/1_LUAD/LUAD_casanovo+pXg")
samples <- c("VY2021_C3N-02145.T", "VY2021_C3N-01416.T", "VY2021_C3N-01024.T",
             "VY2021_C3N-01016.T", "VY2021_C3N-00579.T", "VY2021_C3N-00547.T",
             "VY2021_C3N-00199.T", "VY2021_C3N-00169.T", "VY2021_C3L-02549.T",
             "VY2021_C3L-01632.T", "KZ2022_C3N-02423.T", "KZ2022_C3N-01488.T",
             "KZ2022_C3N-01414.T", "KZ2022_C3N-01030.T", "KZ2022_C3N-01023.T",
             "KZ2022_C3N-00556.T", "KZ2022_C3N-00169.T", "KZ2022_C3L-02549.T",
             "KZ2022_C3L-02348.T", "KZ2022_C3L-00973.T")

data_summary <- as.data.frame(matrix(ncol=7, nrow=0), stringsAsFactors = FALSE, )

for(sample in samples) {
  print(sample)
  
  study <- strsplit(sample, "_")[[1]][1]
  sample_name <- strsplit(sample, "_")[[1]][2]
  
  CSNV <- read.csv(
    paste(sample, ".cnsv.pXg.feat.fdr.ba", sep = ""),
    header = T,
    sep="\t"
  )
  
  data_summary <- rbind(data_summary, 
                        c(sample_name, sample, "Canonical", 
                          get_number_of_psms(CSNV, "true"), get_ba_ratio_of_psms(CSNV, "true"),
                          get_number_of_peptides(CSNV, "true"), get_ba_ratio_of_peptides(CSNV, "true")))
  data_summary <- rbind(data_summary, 
                        c(sample_name, sample, "Noncanonical", 
                          get_number_of_psms(CSNV, "false"), get_ba_ratio_of_psms(CSNV, "false"),
                          get_number_of_peptides(CSNV, "false"), get_ba_ratio_of_peptides(CSNV, "false")))
  
}


colnames(data_summary) <- c("Sample_name", "Sample", "Class", "No. PSMs", "Binder ratio PSMs", "No. Pepts", "Binder ratio Pepts")
colnames(data_summary)

data_summary$`No. PSMs` <- as.numeric(data_summary$`No. PSMs`)
data_summary$`Binder ratio PSMs` <- as.numeric(data_summary$`Binder ratio PSMs`)
data_summary$`No. Pepts` <- as.numeric(data_summary$`No. Pepts`)
data_summary$`Binder ratio Pepts` <- as.numeric(data_summary$`Binder ratio Pepts`)

plot <- ggplot(data_summary[data_summary$Class == "Noncanonical",], 
               aes(x = `Binder ratio PSMs`, y = `No. PSMs`, color = factor(Sample))) +
  geom_point(size = 3) +
  labs(title = "Performance test for non-canonical PSMs") +
  theme_minimal()

ggsave("performance_noncanonical.png", plot, width = 8, height = 6, dpi = 600)

plot <- ggplot(data_summary[data_summary$Class == "Canonical",], 
               aes(x = `Binder ratio Pepts`, y = `No. Pepts`, color = factor(Sample))) +
  geom_point(size = 3) +
  labs(title = "Performance test for canonical peptides") +
  theme_minimal()

ggsave("performance_canonical_peptides.png", plot, width = 8, height = 6, dpi = 600)





