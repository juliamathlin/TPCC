library("sivs")
library("glmnet")
library("pROC")
library("ggplot2")

load("./BCR/2_Ratio_Data_BCR_Groups_DE_ROTS.RData")

file <- data.use.norm

# get genes that have fdr less that 0.05 in bcr groups 
genes <- rots.res.df[which(rots.res.df$fdr < 0.05), ]
df_filtered <- file[rownames(file) %in% rownames(genes), ]


#get bcr status
clinical <- read.csv("./Clinical.csv", sep=";")
rownames(clinical) <- as.character(clinical[,1])
clinica <- clinical[, c("id", "bcr", "Tissue")]
# Remove duplicate rows 
#bcrStatus1 <- unique(bcrStatus1)
rownames(clinica) <- clinica[,"id"] 
bcr <- clinica["bcr"]
Tissue <- clinica["Tissue"]

#transpose df_filtered
df_filtered_t <- t(df_filtered)

merged_df1 <- merge(df_filtered_t, bcr, by = 'row.names', all = TRUE)
rownames(merged_df1) <- merged_df1[,1]
merged_df1 <- subset(merged_df1, select = -Row.names)

merged_df <- merge(merged_df1, Tissue, by = 'row.names', all = TRUE)
rownames(merged_df) <- merged_df[,1]
merged_df <- subset(merged_df, select = -Row.names)


cancer <- merged_df[merged_df[ ,"Tissue"] == "Cancer", ]
DiseaseOutcomeVector <- cancer$bcr
cancer <- subset(cancer, select = -c(bcr, Tissue))
benign <- merged_df[merged_df[ ,"Tissue"] == "Benign", ]
benign <- subset(benign, select = -c(bcr, Tissue))

data.benSivs=log2(benign+1) 
data.canSivs=log2(cancer+1) 
cancerMSivs <- as.matrix(data.canSivs)
benignMSivs <- as.matrix(data.benSivs)
new_dfSivs <- cancerMSivs-benignMSivs

# remove columns that have NA
merged_df <- merged_df[complete.cases(merged_df), ]

sivs_obj <- sivs::sivs(x = new_dfSivs, y = factor(DiseaseOutcomeVector), verbose = "detailed", progressbar = TRUE, family = "binomial")

# Plot
layout(mat = matrix(c(1,2),
                    nrow = 1,
                    byrow = T))
{
  plot(sivs_obj)
  layout(1)
}
features <- sivs::suggest(sivs_obj, strictness =0.5)
features
plot(sivs_obj)
