# AUC estimates and model selection apoptosis index

library(knitr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(glmnet)
library(plotly)
library(gplots)
library(viridis)

colours <- c("#EA7580","#D8C99E","#172869","#019875")

source("../help_plot_AUC2.R")

#Prep data on model outcomes for plotting
res_el <- readRDS("./12-output/1-output/res_el.Rds")
el1 <- combine_sets(res_el, "Elastic Net - 1", type = "slope")

res_el <- readRDS("./12-output/2-output/res_el.Rds")
el2 <- combine_sets(res_el, "Elastic Net - 2", type = "slope")

res_el <- readRDS("./12-output/3-output/res_el.Rds")
el3 <- combine_sets(res_el, "Elastic Net - 3", type = "slope")

res_el <- readRDS("./12-output/4-output/res_el.Rds")
el4 <- combine_sets(res_el, "Elastic Net - 4", type = "slope")

pdat <- rbind(el1,el2,el3,el4)
pdat <- as.data.frame(pdat)

pdat$n <- as.numeric(pdat$n)
pdat$AUC_val <- as.numeric(pdat$AUC_val)
pdat$mean_deviance <- as.numeric(pdat$AUC_oob) # sample size to small for type.measure=AUC
pdat$AUC_tr <- as.numeric(pdat$AUC_tr)
pdat$Slope <- as.numeric(pdat$Slope)
pdat$Intercept <- as.numeric(pdat$Intercept)
pdat$Type <- factor(pdat$Type)

d1 <- pdat %>%
  ggplot() +
  geom_line(aes(x = n,
                y = mean_deviance,
                colour = Type)) +
  xlab("Number of input CpGs") +
  ylab("mean deviance (across Lambda)") +
  labs(color="")+
  theme_minimal() +
  scale_colour_manual(values = colours,
                      aesthetics = "colour")

d1 <- ggplotly(d1)

d2 <- pdat %>%
  ggplot(aes(x = Slope,
             y = mean_deviance,
             colour = Type,
             text = n)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  scale_colour_manual(values = colours,
                      aesthetics = "colour")+
  xlab("Unitless slope calibration plot (Cox)") +
  ylab("Mean deviance") +
  labs(color="") +
  xlim(-0.1,1.5) 

d2 <- ggplotly(d2)

d3 <- pdat %>%
  ggplot(aes(x = Intercept,
             y = Slope,
             colour = Type,
             text = n)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  scale_colour_manual(values = colours,
                      aesthetics = "colour")+
  xlab("Calibration intercept") +
  ylab("Unitless slope calibration plot (Cox)") +
  labs(color="") +
  xlim(-10,10) +
  ylim(-0.1,1.5)

d3 <- ggplotly(d3)

d4 <- pdat %>%
  ggplot(aes(x = Intercept,
             y = mean_deviance,
             colour = Type,
             text = n)) +
  geom_point(size = 0.5) +
  theme_minimal() +
  scale_colour_manual(values = colours,
                      aesthetics = "colour")+
  xlab("Calibration intercept") +
  ylab("Mean deviance") +
  xlim(-10,10) + 
  labs(color="") 

d4 <- ggplotly(d4)

subplot(style(d1, showlegend=T), 
        style(d2, showlegend=F),
        style(d3, showlegend=F),
        style(d4, showlegend=F), 
        nrows = 2, titleX=T, titleY=T, margin = 0.07)

# Selecting Approach 2: 3000 Input Cpgs, after penalizing Elastic net, gives 114 CpGs
# Not in plot bounderies with slope/intercept: Elastic Net -3, n=10: AUC_val=1, Slope = 297, Intercept = 118, Mean deviance = 1.23


