library(pheatmap)

load("./BCR/2_Ratio_Data_BCR_Groups_DE_ROTS.RData")

# use normalized data
df_hm <- data.norm
x <- as.matrix(df_hm[, 1:ncol(df_hm)])
x <- apply(x, 2, as.numeric)  
rownames(x) <- rownames(df_hm)

# Clinical data
clin <- clinical
clin$id <- NULL
clin$Margin <- NULL
clin$PSA <- NULL
clin$PSA_original <- NULL

clin$bcr <- as.factor(clin$bcr)
clin$pN <- as.character(clin$pN)
clin$GGG <- as.factor(clin$GGG)
clin$Tissue <- as.factor(clin$Tissue)

# Malignant vs Benign Heatmap

annotation_colors <- list(
  bcr = c("#FF964D", "#57B957"),    
  pN = c("1" = "deeppink", "2"="deepskyblue", "NA"="gray"), 
  GGG =  colorRampPalette(c("#FFF1E6", "#A04000"))(5),
  Tissue = c("darkblue", "#CF9FFF")
)

levels_bcr <- unique(clin$bcr) # Assuming bcr has distinct values like 0, 1
names(annotation_colors$bcr) <- levels_bcr

levels_pN <- unique(clin$pN) # Assuming pn has distinct values like 1, 2, 3
names(annotation_colors$pN) <- c("1", "0", "x")

levels_GGG <- unique(clin$GGG) # Assuming pn has distinct values like 1, 2, 3
names(annotation_colors$GGG) <- c(1, 2, 3, 4, 5)

levels_Tissue <- unique(clin$Tissue) # Assuming pn has distinct values like 1, 2, 3
names(annotation_colors$Tissue) <- c("Cancer", "Benign")

scaleData<-function(x){
  require(matrixStats)
  rmean <- rowMeans(x,na.rm = TRUE)
  rsd <- rowSds(x,na.rm = TRUE)
  xS <- (x-rmean)/rsd
  return(xS)
}

drawHeatmap <- function(x.s, x.t , group,  gcol, clustMethod, bk2,  filename, fileType, w = 5,  he = 10 ){
  
  require(gplots)
  require(RColorBrewer)
  require(pheatmap)
  hcdR <- hclust(dist(as.matrix(x.s), 'euclidean'), method =  clustMethod) 
  hcdC <- hclust(dist(as.matrix(t(x.s)), 'euclidean'), method =  clustMethod)
  cols = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(500)
  if(fileType == 'pdf'){
    pdf(paste0(filename, '.', fileType), width = w, height = he)
  }
  if(fileType == 'tiff'){
    tiff(paste0(filename, '.', fileType), width = w, height = he, units = 'in', res = 3000)
  }
  if(fileType == 'svg'){
    svg(paste0(filename, '.', fileType), width = w, height = he)
  }
  #
  pheatmap(x.t, scale = 'none', cluster_cols = hcdC, 
           cluster_rows = hcdR, annotation_colors = annotation_colors,
           main = '', fontsize = 6, 
           border_color = NA, 
           fontsize_row = 0.01, 
           #fontsize_col = 4,
           color=cols,
           gaps_col = 2,
           breaks = bk2,
           cellwidth = 3, cellheight = 0.02, 
           treeheight_row = 0, 
           annotation_col = clin, 
           filename = NA)
  dev.off()
}


data.scaled <- scaleData(x)
x2 <- data.scaled
x2[x2 > 2] <- 2
bk2 = unique(c(seq(-2, 2, length = 500)))


drawHeatmap(data.scaled, x2, group, gcol=annotation_colors, 'ward.D', bk2, 'HeatmapBvsM', 'pdf')

####################
# BCR heatmap
####################

# use normalized data
file <- data.log.ratio

# get genes that have fdr less that 0.05 in bcr groups 
genes <- rots.res.df[which(rots.res.df$fdr < 0.05 & (rots.res.df$logfc > 1 | rots.res.df$logfc < -1)), ]
#get only columns that are of the significant genes of bcr groups - make sure there are no NA or rows full of zero
df_filtered <- file[rownames(file) %in% rownames(genes), ]
rownames(df_filtered) <- gsub("\\s*\\*\\* LowExpr\\s*", "", rownames(df_filtered))

x <- data.matrix(df_filtered[,1:ncol(df_filtered)])
rnames <- rownames(df_filtered)                        
rownames(x) <- rnames                 

# Add clinical data
clin <- clinical
clin$id <- NULL
clin$Tissue <- NULL
clin$PSA_original <- NULL

clin$bcr <- as.factor(clin$bcr)
clin$pN <- as.character(clin$pN)
clin$GGG <- as.factor(clin$GGG)
clin$Margin <- as.factor(clin$Margin)
clin$PSA <- factor(clin$PSA)


annotation_colors_bcr <- list(
  bcr = c("#FF964D", "#57B957"),      
  pN = c("1" = "deeppink", "2"="deepskyblue", "NA"="gray"), 
  GGG =  colorRampPalette(c("#FFF1E6", "#A04000"))(5),
  PSA = colorRampPalette(c("honeydew", "darkgreen"))(6),
  Margin = c("purple", "#40E0D0"))

levels_bcr <- unique(clin$bcr) # Assuming bcr has distinct values like 0, 1
names(annotation_colors_bcr$bcr) <- levels_bcr

levels_pN <- unique(clin$pN) # Assuming pn has distinct values like 1, 2, 3
names(annotation_colors_bcr$pN) <- c("1", "0", "x")

levels_GGG <- unique(clin$GGG) # Assuming pn has distinct values like 1, 2, 3
names(annotation_colors_bcr$GGG) <- c(1, 2, 3, 4, 5)

levels_PSA <- unique(clin$PSA) 
names(annotation_colors_bcr$PSA) <- c("<5", "<10", "<15", "<20", "<30",">30")

levels_Margin <- unique(clin$Margin) 
names(annotation_colors_bcr$Margin) <- c("1", "0")


clin <- clin[rownames(clin) %in% colnames(df_filtered), ]
clin <- clin[match(colnames(df_filtered), rownames(clin)), ]

# make a gab according to clusters
ord <- order(clin$bcr)
df_filtered_ord <- df_filtered[, ord]
clinical_ord <- clin[ord, , drop = FALSE]
gaps_col <- which(diff(as.numeric(clinical_ord$bcr)) != 0)

mx <- max(abs(range(df_filtered, na.rm = TRUE)))
bk <- seq(-mx, mx, length.out = 501)

# Reverse RdBu so NEGATIVE = BLUE, POSITIVE = RED
cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))(500))

svg("Heatmap-BCR.svg", width = 12, height = 14)

pheatmap(
  df_filtered_ord,
  scale = "none",
  cluster_cols = hclust(dist(t(df_filtered_ord), "euclidean"), method = "ward.D"),
  cluster_rows = hclust(dist(df_filtered_ord, "euclidean"), method = "ward.D"),
  annotation_col = clinical_ord,
  annotation_colors = annotation_colors_bcr,
  cutree_cols = 2,  
  color  = cols,
  breaks = bk,
  border_color = NA,
  show_colnames = FALSE,
  treeheight_row = 0,
  treeheight_col = 40,
  cellwidth = 10,
  cellheight = 8
)

dev.off()
