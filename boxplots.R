# Read the file containing genes and FDR values
library(readxl)
library(stringr)
library(ggplot2)

load("./BCR/2_Ratio_Data_BCR_Groups_DE_ROTS.RData")

data <- t(as.data.frame(cpm.data))
data[is.na(data)] <- 0

# Filter genes based on FDR threshold
gene_fdr_data <- data.log.ratio
colnames(gene_fdr_data)[1] <- "Gene"
significant_genes <- gene_fdr_data$Gene[gene_fdr_data$fdr < 0.05]
significant_genes <- gsub(" \\*\\* LowExpr", "", significant_genes)


for (i in significant_genes) {
  benign <- data.frame(data[clinical$Tissue == "Benign", i, drop = FALSE])
  benign[,1] <- log2(benign[,1])
  benignt<- cbind(benign, sampleType="Benign")
  
  malign <- data.frame(data[clinical$Tissue == "Cancer", i, drop = FALSE])
  malign[,1] <- log2(malign[,1])
  malignt<- cbind(malign, sampleType="Malignant")
  
  bindded <- rbind(benignt, malignt)
  
  bindded$paired = row.names(bindded)  
  bindded$paired <- sub("_.*$", "", bindded$paired)
  
  bindded$bcr <- clinical$bcr[
    match(
      bindded$paired,
      sub("_.*$", "", rownames(clinical))
    )
  ]
  
  bindded$bcr <- ifelse(bindded$bcr == 1,
                        "Yes_bcr",
                        "No_bcr")
  # Create ggplot
  p <- ggplot(bindded, aes_string(x = "sampleType", y = i, fill = "sampleType")) +
    geom_boxplot(lwd = 0.8, color = "black") +
    scale_fill_manual(values = c(
      "Malignant" = "#D0E4FF",  # light gray
      "Benign"  = "#F5F5DC"   # beige
    )) +
    geom_point() +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90"),           
      panel.grid.minor = element_line(color = "gray95"),
      axis.text = element_text(size = 40),  
    ) +
    scale_y_continuous(limits = c(-1, 8)) + 
    geom_line(aes(group = paired, color = bcr), lwd = 1) +
    scale_color_manual(values = c("Yes_bcr" = "#FF964D",  
               "No_bcr" = "#57B957")) 
 
  ggsave(paste0("Boxplot_bcr",i,".pdf"), plot = p, width = 10, height = 12)
}

