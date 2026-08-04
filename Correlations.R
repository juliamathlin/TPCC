library(igraph)
library(ggraph)
library(ggrepel)

load("./BCR/2_Ratio_Data_BCR_Groups_DE_ROTS.RData")


datacor <- data.use.norm
rownames(datacor) <- gsub("\\s*\\*\\* LowExpr\\s*", "", rownames(datacor))

geneset1 <- c("BPIFB2", "C11orf96", "NR4A2", "DUSP5", "NR4A3")

geneset2 <- rots.res.df[rots.res.df$fdr < 0.05, ]
geneset2 <- geneset2[!rownames(geneset2) %in% c("CSRNP1", "PAQR8"), ]
rownames(geneset2) <- gsub("\\s*\\*\\* LowExpr\\s*", "", rownames(geneset2))
geneset2 <- rownames(geneset2)

# loop though BCRs and Tissues

conditions <- expand.grid(
  BCR = c(1, 0),
  Tissue = c("Cancer", "Benign"),
  stringsAsFactors = FALSE
)

for(i in seq_len(nrow(conditions))) {
  
  bcr_group <- conditions$BCR[i]
  tissue_group <- conditions$Tissue[i]
  
  analysis_name <- paste0(
    "BCR",
    bcr_group,
    "_",
    toupper(tissue_group)
  )
  
samples_keep <- rownames(clinical)[
    clinical$bcr == bcr_group &
      clinical$Tissue == tissue_group
  ]
  
  df <- datacor[, colnames(datacor) %in% samples_keep, drop = FALSE]


# correlate  
  cor_results <- data.frame()
  
  for(gene_x in geneset1){
    
    if(!gene_x %in% rownames(df))
      next
    
    for(gene_y in geneset2){
      
      if(!gene_y %in% rownames(df))
        next
      
      datacor_x <- as.numeric(df[gene_x, ])
      datacor_y <- as.numeric(df[gene_y, ])
      
      valid <- is.finite(datacor_x) & is.finite(datacor_y)
      
      if(sum(valid) < 3)
        next
      
      ct <- tryCatch(
        cor.test(datacor_x[valid], datacor_y[valid], method = "pearson"),
        error = function(e) NULL
      )
      
      if(!is.null(ct)){
        
        cor_results <- rbind(
          cor_results,
          data.frame(
            Gene1 = gene_x,
            Gene2 = gene_y,
            Correlation = unname(ct$estimate),
            PValue = ct$p.value
          )
        )
      }
    }
  }
  
  if(nrow(cor_results) == 0){
    message("No correlations found for ", analysis_name)
    next
  }
  
  
  write.csv(
    cor_results,
    paste0("correlation_results_", analysis_name, ".csv"),
    row.names = FALSE
  )
  
# graphs
  df_filtered <- subset(
    cor_results,
    abs(Correlation) > 0.6 &
      PValue < 0.05
  )
  
  df_filtered$EdgeColor <- ifelse(
    df_filtered$Correlation > 0,
    "deepskyblue2",
    "brown"
  )
  
  graph <- graph_from_data_frame(
    df_filtered,
    directed = FALSE
  )
  
  layout_positions <- create_layout(
    graph,
    layout = "fr"
  )
  
   layout_positions$highlight <- ifelse(
    layout_positions$name %in% geneset1,
    "darkgreen",
    "darkgray"
  )
  
  p <- ggraph(layout_positions) +
    geom_edge_link(
      aes(color = EdgeColor),
      linewidth = 1,
      edge_alpha = 1
    ) +
    geom_node_point(
      aes(color = highlight),
      size = 10
    ) +
    geom_text_repel(
      aes(x = x, y = y, label = name),
      data = layout_positions,
      size = 20,
      seed = 123
    ) +
    scale_edge_color_identity() +
    scale_color_identity() +
    theme_void() +
    ggtitle(analysis_name)
  
  ggsave(
    paste0("network_", analysis_name, ".pdf"),
    plot = p,
    width = 20,
    height = 20
  )
}
