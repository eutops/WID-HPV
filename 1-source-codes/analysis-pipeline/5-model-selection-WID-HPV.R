# AUC estimates and model selection WID-HPV

library(knitr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(glmnet)
library(plotly)
library(gplots)
library(viridis)

colour.pal.d8 <- c("#EA7580","#F2949A","#F6B3A1","#D8C99E","#3AB9AC","#109DB7","#0C6EA5","#172869")
colour.pal.c <- colorRampPalette(colour.pal.d8)

source("../help_plot_AUC.R")

## OOB AUC

load("./4-output/res_ridge_f.Rdata")
load("./4-output/res_el.Rdata")
load("./4-output/res_lasso_f.Rdata")

ridge <- combine_sets(res_ridge, "Ridge", type = "slope")
lasso <- combine_sets(res_lasso, "Lasso", type = "slope")
el <- combine_sets(res_el, "Elastic Net", type = "slope")

res_ridge_epi <- readRDS("./4-output/res_ridge_epi.Rds")
res_lasso_epi <- readRDS("./4-output/res_lasso_epi.Rds")
res_ridge_imm <- readRDS("./4-output/res_ridge_imm.Rds")
res_lasso_imm <- readRDS("./4-output/res_lasso_imm.Rds")
res_ridge_epi_imm <- readRDS("./4-output/res_ridge_epi_imm.Rds")
res_lasso_epi_imm <- readRDS("./4-output/res_lasso_epi_imm.Rds")

ridge_epi <- combine_sets(res_ridge_epi, "Ridge - Epithelial", type = "slope")
lasso_epi <- combine_sets(res_lasso_epi, "Lasso - Epithelial", type = "slope")
ridge_imm <- combine_sets(res_ridge_imm, "Ridge - Immune", type = "slope")
lasso_imm <- combine_sets(res_lasso_imm, "Lasso - Immune", type = "slope")
ridge_comb <- combine_sets(res_ridge_epi_imm, "Ridge - Combined", type = "slope")
lasso_comb <- combine_sets(res_lasso_epi_imm, "Lasso - Combined", type = "slope")

pdat <- rbind(ridge, lasso, el, ridge_epi, lasso_epi,  ridge_imm, lasso_imm, ridge_comb, lasso_comb)

pdat <- as.data.frame(pdat)

pdat$n <- as.numeric(pdat$n)
pdat$AUC_val <- as.numeric(pdat$AUC_val)
pdat$AUC_oob <- as.numeric(pdat$AUC_oob)
pdat$AUC_tr <- as.numeric(pdat$AUC_tr)
pdat$Slope <- as.numeric(pdat$Slope)
pdat$Type <- factor(pdat$Type, levels = c("Ridge", "Lasso", "Elastic Net", "Ridge - Epithelial", "Lasso - Epithelial", "Ridge - Immune", "Lasso - Immune", "Ridge - Combined", "Lasso - Combined"))

ggplotly(plot_performance_compare(pdat, pdat$AUC_oob, sets = 12))


## Slope

d <- pdat %>%
  ggplot(aes(x = Slope,
             y = AUC_oob,
             colour = Type,
             text = n)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  xlab("Slope") +
  ylab("Out-of-bag AUC") +
  xlim(-0.1,1.5)

ggplotly(d)

#Selecting lasso - epithelial n = 500; returns n = 27 CpGs