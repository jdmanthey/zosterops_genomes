# Last editted by Ethan Gyllenhaal, 2 December 2024
# Code for re-analyzing Zosterops genetic diversity data
# Also re-made plots separating by island group
# Some edits made to text in Illustrator

library(ggplot2)
library(gridExtra)
library(cowplot)
library(tidyverse)
library(glmnet)
library(corrplot)

setwd('/path/to/directory')

data <- read.csv("zosterops_genomes_appendix_2024.csv") %>%
  filter(Island!="Not in Solomons")

data_nogizo <- data %>%
  filter(Island!="Gizo")

colors <- c("#059948", "#7B499D", "#b3a68f", "#F26522")
shapes <- c(16, 17, 15, 18)

# Figure 2 plots ####

# island size vs recent Ne (legend adjusted in illustrator)
i_neR <- ggplot(data=data, aes(x=Island.Size..log., y=MSMC.Recent.Pop.Size)) +
  geom_smooth(method="lm", se=F, show.legend = FALSE, color="black") + 
  geom_point(size=2.5, show.legend = TRUE, aes(color=Island.Group, shape=Island.Group)) +
  scale_color_manual(values = colors) +  scale_shape_manual(values = shapes) + theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle=90, hjust=0.5),
        legend.position = "inside", legend.position.inside = c(0.25, 0.75), legend.key.size = unit(1,'mm'), legend.title = element_text(size=9))

# island size vs historic Ne
i_neH <- ggplot(data=data, aes(x=Island.Size..log., y=MSMC.Harmonic.Mean.Pop.Size)) +
  geom_smooth(method="lm", se=F, show.legend = FALSE, color="black") +
  geom_point(size=2.5, show.legend = FALSE, aes(color=Island.Group, shape=Island.Group)) +
  scale_color_manual(values = colors) +  scale_shape_manual(values = shapes) + theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle=90, hjust=0.5))

# island size vs heterozygosity
i_div <- ggplot(data=data, aes(x=Island.Size..log., y=obs_het)) +
  geom_smooth(method="lm", se=F, show.legend = FALSE, color="black") + ylim(0,0.0065) +
  geom_point(size=2.5, show.legend=FALSE, aes(color=Island.Group, shape=Island.Group)) +
  scale_color_manual(values = colors) +  scale_shape_manual(values = shapes) + theme_bw() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.text.y = element_text(angle=90, hjust=0.5))

# island size vs heterozygosity SD
i_sd <- ggplot(data=data, aes(x=Island.Size..log., y=het_sd)) +
  geom_smooth(method="lm", se=F, show.legend = FALSE, color="black") +
  geom_point(size=2.5, show.legend = FALSE, aes(color=Island.Group, shape=Island.Group)) +
  scale_color_manual(values = colors) +  scale_shape_manual(values = shapes) +  theme_bw() +
  theme(axis.text.y = element_text(angle=90, hjust=0.5))

# island size vs ERV fixation
i_erv <- ggplot(data=data, aes(x=Island.Size..log., y=NR.Homozygous.TEs)) +
  geom_smooth(method="lm", se=F, show.legend = FALSE, color="black") +
  geom_point(size=2.5, show.legend = FALSE, aes(color=Island.Group, shape=Island.Group)) +
  scale_color_manual(values = colors) +  scale_shape_manual(values = shapes) +  theme_bw() +
  theme(axis.text.y = element_text(angle=90, hjust=0.5))

# island size vs ROH
i_roh <- ggplot(data=data, aes(x=Island.Size..log., y=Length.ROH/1000000)) +
  geom_smooth(method="lm", se=F, show.legend = FALSE, color="black") +
  geom_point(size=2.5, show.legend = FALSE, aes(color=Island.Group, shape=Island.Group)) +
  scale_color_manual(values = colors) +  scale_shape_manual(values = shapes) +  theme_bw() +
  theme(axis.text.y = element_text(angle=90, hjust=0.5))

# plot together
plot_grid(i_neR, i_neH, i_div, i_sd, i_erv, i_roh, nrow=2, 
             rel_heights = c(0.45, 0.55))

# simple linear models for individual Rsq
summary(lm(MSMC.Recent.Pop.Size~Island.Size..log., data=data))
summary(lm(MSMC.Harmonic.Mean.Pop.Size~Island.Size..log., data=data))
summary(lm(obs_het~Island.Size..log., data=data))
summary(lm(het_sd~Island.Size..log., data=data))
summary(lm(NR.Homozygous.TEs~Island.Size..log., data=data))
summary(lm(Length.ROH~Island.Size..log., data=data))
summary(lm(Length.ROH~Island.Size..log., data=data_nogizo))


# Correlogram ####

# make matrices for input
subset_data <- data.matrix(data[,c('Island.Size..log.', 'MSMC.Harmonic.Mean.Pop.Size', 'MSMC.Recent.Pop.Size', 'obs_het', 'het_sd', 'NR.Homozygous.TEs', 'Length.ROH', 'X..ROH')])
cor_matrix = cor(subset_data)
subset_nogizo <- data.matrix(data_nogizo[,c('Island.Size..log.', 'MSMC.Harmonic.Mean.Pop.Size', 'MSMC.Recent.Pop.Size', 'obs_het', 'het_sd', 'NR.Homozygous.TEs', 'Length.ROH', 'X..ROH')])
cor_nogizo = cor(subset_nogizo)

# main one for labels
corrplot(cor_matrix, method="square", addCoef.col = "black", tl.col = "black", tl.srt = 45)
# plot others for smaller text coefficients
# note that this leads to some stray stuff from the original one to clean up
corrplot(cor_matrix, type="lower", method="square", addCoef.col = "black", tl.pos = "n", add=T, number.cex=0.8)
corrplot(cor_nogizo, type="upper", method="square", addCoef.col = "black", tl.pos = "n", add=T, number.cex=0.8)

# ROH plot ####
# calculate ratio in 2 ways, same result, plot the one with Gizo
mean_ratio <- mean(data$Length.ROH/data$X..ROH)
mean_ratio_nogizo <- mean(data_nogizo$Length.ROH/data_nogizo$X..ROH)
# make a plot with a line for the mean ratio
roh_plot <- ggplot(data=data, aes(x=Length.ROH, y=X..ROH)) +
  geom_point(size=2.5, show.legend = TRUE, aes(color=Island.Group, shape=Island.Group)) +
  scale_color_manual(values = colors) +  scale_shape_manual(values = shapes) +  theme_bw() +
  theme(axis.text.y = element_text(angle=90, hjust=0.5), legend.position = "inside", legend.position.inside = c(0.15, 0.85), legend.key.size = unit(5,'mm'), legend.title = element_text(size=10)) +
  geom_abline(slope=1/mean_ratio, intercept=0)
roh_plot

# LASSO ####
# note that Rsq was calculated, but interpretation is controvesial, so left it out
# (they are all high)

# TEs
y_TE <- data$NR.Homozygous.TEs
full_x_TE <- data.matrix(data[,c('Island.Size..log.', 'obs_het', 'het_sd', 'Length.ROH', 'MSMC.Harmonic.Mean.Pop.Size', 'MSMC.Recent.Pop.Size', 'X..ROH')])
cv_full_TE <- cv.glmnet(x=full_x_TE, y=y_TE, alpha=1, grouped=FALSE)
plot(cv_full_TE)
best_lambda_full_TE <- cv_full_TE$lambda.min
best_lambda_full_TE
best_full_TE <- glmnet(x=full_x_TE, y=y_TE, alpha=1, lambda = best_lambda_full_TE)
coef(best_full_TE)
pred_full_TE <- predict(best_full_TE, s=best_lambda_full_TE, newx = full_x_TE)
sst_full_TE <- sum((y_TE - mean(y_TE))^2)
sse_full_TE <- sum((pred_full_TE - y_TE)^2)
rsq_full_TE <- 1 - sse_full_TE/sst_full_TE
rsq_full_TE

# island size
y_island <- data$Island.Size..log.
full_x_island <- data.matrix(data[,c('NR.Homozygous.TEs', 'obs_het', 'het_sd', 'Length.ROH', 'MSMC.Harmonic.Mean.Pop.Size', 'MSMC.Recent.Pop.Size', 'X..ROH')])
cv_full_island <- cv.glmnet(x=full_x_island, y=y_island, alpha=1, grouped=FALSE)
plot(cv_full_island)
best_lambda_full_island <- cv_full_island$lambda.min
best_lambda_full_island
best_full_island <- glmnet(x=full_x_island, y=y_island, alpha=1, lambda = best_lambda_full_island)
coef(best_full_island)
pred_full_island <- predict(best_full_island, s=best_lambda_full_island, newx = full_x_island)
sst_full_island <- sum((y_island - mean(y_island))^2)
sse_full_island <- sum((pred_full_island - y_island)^2)
rsq_full_island <- 1 - sse_full_island/sst_full_island
rsq_full_island

## not used

# island size without gizo
y_island_nogizo <- data$Island.Size..log.
full_x_island_nogizo <- data.matrix(data[,c('NR.Homozygous.TEs', 'obs_het', 'het_sd', 'Length.ROH', 'MSMC.Harmonic.Mean.Pop.Size', 'MSMC.Recent.Pop.Size', 'X..ROH')])
cv_full_island_nogizo <- cv.glmnet(x=full_x_island_nogizo, y=y_island_nogizo, alpha=1, grouped=FALSE)
plot(cv_full_island_nogizo)
best_lambda_full_island_nogizo <- cv_full_island_nogizo$lambda.min
best_lambda_full_island_nogizo
best_full_island_nogizo <- glmnet(x=full_x_island_nogizo, y=y_island_nogizo, alpha=1, lambda = best_lambda_full_island_nogizo)
coef(best_full_island_nogizo)
pred_full_island_nogizo <- predict(best_full_island_nogizo, s=best_lambda_full_island_nogizo, newx = full_x_island_nogizo)
sst_full_island_nogizo <- sum((y_island_nogizo - mean(y_island_nogizo))^2)
sse_full_island_nogizo <- sum((pred_full_island_nogizo - y_island_nogizo)^2)
rsq_full_island_nogizo <- 1 - sse_full_island_nogizo/sst_full_island_nogizo
rsq_full_island_nogizo


# heterozygosity
y_het <- data$obs_het
full_x_het <- data.matrix(data[,c('NR.Homozygous.TEs', 'Island.Size..log.', 'het_sd', 'Length.ROH', 'MSMC.Harmonic.Mean.Pop.Size', 'MSMC.Recent.Pop.Size')])
cv_full_het <- cv.glmnet(x=full_x_het, y=y_het, alpha=1, grouped=FALSE)
plot(cv_full_het)
best_lambda_full_het <- cv_full_het$lambda.min
best_lambda_full_het
best_full_het <- glmnet(x=full_x_het, y=y_het, alpha=1, lambda = best_lambda_full_het)
coef(best_full_het)
pred_full_het <- predict(best_full_het, s=best_lambda_full_het, newx = full_x_het)
sst_full_het <- sum((y_het - mean(y_het))^2)
sse_full_het <- sum((pred_full_het - y_het)^2)
rsq_full_het <- 1 - sse_full_het/sst_full_het
rsq_full_het


# TE vs Div, not actually used

te_div <- ggplot(data=data, aes(x=obs_het, y=NR.Homozygous.TEs)) +
  geom_smooth(method="lm", se=T, show.legend = FALSE, color="black") +
  geom_point(size=2, show.legend = FALSE, aes(color=Island.Group, shape=Island.Group)) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  theme_bw()
te_div
