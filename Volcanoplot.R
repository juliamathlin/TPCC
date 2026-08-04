library(ggplot2)
library(ggrepel)
library(dplyr)
library(svglite)
library(EnhancedVolcano)

load("BvsM_Paired_Analyses/1_PreProcessing_Done_ROTS_DE_Done.RData")

# select file
volcanoData <- rots.res.df %>%
  select(fdr,
         logfc,
         pvalue) 

volcanoData <- volcanoData[complete.cases(volcanoData),]

# select genes
genelist <- c("MYC", "TP53", "AR", "SMAD4", "ERG")


svglite(file = "VolcanoPlot_BvsM_allGenes.svg",
        width = 7,
        height = 10
) 

rearranging <- as.character(rownames(volcanoData)) %in% genelist
volcanoData <- rbind(volcanoData[!rearranging,], volcanoData[rearranging,])


keyvalsHvsL <- ifelse(volcanoData$logfc < -log2(1.5) & volcanoData$fdr < 0.01, '#92c5de',
                      ifelse(volcanoData$logfc > log2(1.5) & volcanoData$fdr < 0.01, '#f4a582',
                             ifelse(rownames(volcanoData) %in% "keyvalsHvsL", "red",
                                    'grey')
                      )
)

keyvalsHvsL[volcanoData$symbol %in% genelist] <- 'red'
keyvalsHvsL[is.na(keyvalsHvsL)] <- 'grey'
names(keyvalsHvsL)[keyvalsHvsL == '#f4a582'] <- 'high'
names(keyvalsHvsL)[keyvalsHvsL == 'grey'] <- 'mid'
names(keyvalsHvsL)[keyvalsHvsL == '#92c5de'] <- 'low'
names(keyvalsHvsL)[keyvalsHvsL == 'red'] <- 'Genes'
names(keyvalsHvsL)[keyvalsHvsL == 'grey'] <- 'NA'



EnhancedVolcano(volcanoData,
                lab = rownames(volcanoData),
                x = "logfc",
                y = "pvalue",
                selectLab = genelist,
                # pCutoff = 10e-14,
                #FCcutoff = 1.5,
                pointSize = 2,
                labSize = 6.5,
                labCol = 'black',
                labFace = 'italic',
                boxedLabels = FALSE,
                colAlpha = 4/5,
                legendPosition = 'none',
                #legendLabSize = 5,
                #legendIconSize = 3,
                endsConnectors="last",
                drawConnectors = TRUE,
                typeConnectors="open",
                widthConnectors = 0.001,
                #max.overlaps = "Inf",
                colConnectors = 'red',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                cutoffLineType = 'blank',
                colCustom = keyvalsHvsL,
                xlim= c(-6, 6),
                ylim= c(0, 9),
                title = "Benign_vs_Malign_Gene_Expression")
dev.off()

