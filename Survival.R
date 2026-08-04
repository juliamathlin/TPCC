library(survminer)
library(survival)
library(ggplot2)

load("./BCR/2_Ratio_Data_BCR_Groups_DE_ROTS.RData")

data <- as.data.frame(data.log.ratio)
rownames(data) <- gsub(" \\*\\* LowExpr", "", rownames(data))

genes <- c("NR4A3","NR4A2","DUSP5","C11orf96", "BPIFB2")
dataGenes <- t(data[genes, ])
colnames(dataGenes) <- genes

clinic <- cbind(clinical[rownames(dataGenes), , drop = FALSE], dataGenes)

clinic$RALP_date = as.Date(clinic$RALP_date, format = "%d.%m.%Y")
clinic$bcr_date = as.Date(clinic$bcr_date, format = "%d.%m.%Y")
clinic$survival_time = clinic$bcr_date - clinic$RALP_date

# Multivariate Cox for top 5 sivs genes
cox.model <- coxph(Surv(survival_time/365, bcr == "1") ~ NR4A3 + NR4A2 + DUSP5 + C11orf96 + BPIFB2,
                   data = clinic)
summary(cox.model)

# plot
ggforest(cox.model, data = clinic,
         main = "Cox model for log2fc2 genes",
         cpositions = c(0.001, 0.15, 0.5),
         fontsize = 1, 
        # annotate = FALSE)
)
ggsave(filename = "Cox_forest_sivs_genes.pdf", height = 5, width = 5)

##################
# Survival 

# Cox risk score
clinic$riskScore <- predict(cox.model, type = "lp")

clinic$group5 <- ifelse(clinic$riskScore > 0, "High Risk", "Low Risk")
clinic$group5 <- factor(clinic$group5, levels = c("High Risk", "Low Risk"))

med_cut <- median(clinic$riskScore, na.rm = TRUE)
clinic$group_median <- ifelse(clinic$riskScore > med_cut, "High Risk", "Low Risk")

# KM
fitKM5 <- survfit(Surv(survival_time/365, bcr == "1") ~ group5, data = clinic)

# p-value by log-rank test
logrank_test <- survdiff(Surv(survival_time/365, bcr == "1") ~ group5, data = clinic)
p_value <- 1 - pchisq(logrank_test$chisq, df = length(logrank_test$n) - 1)
formatted_pval <- paste0("p = ", signif(p_value, digits = 5))


ggsurvplot(
  fitKM5,
  data = clinic,
  xlab = "Time (years)",
  break.x.by = 1,
  risk.table = TRUE,
  cumevents = TRUE,
  censor = TRUE,
  surv.scale = "percent",
  palette = c("seagreen", "salmon"),
  legend.labs = levels(clinic$group5),
  pval = formatted_pval,
  pval.method = TRUE,
  fontsize = 9,  # affects risk table text slightly
  title = "SIVS + NR4A3",
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE, 
  risk.table.y.text = TRUE,      
  ggtheme = theme_classic(base_size = 30) +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 28, face = "bold"),
      axis.text = element_text(size = 30),
      legend.title = element_text(size = 30, face = "bold"),
      legend.text = element_text(size = 30),
      plot.margin = margin(20, 20, 20, 20)
    )
)


ggsave(filename = "survival.pdf", height = 5, width = 5)
