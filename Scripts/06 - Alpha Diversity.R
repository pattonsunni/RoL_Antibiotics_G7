# Author: Sunni Patton
# Last edited: 06/26/24
# Title: Alpha diversity 

## Set seed ====
set.seed(123)

## Load libraries ====
library(picante)
library(phyloseq)
library(car)
library(stats)
library(ggpubr)
library(DescTools)
library(rstatix)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(pbkrtest)
library(microViz)


## Load data ====
readRDS(here::here("Data/05 - Phyloseq - Output/ps.rare.rds")) -> ps.rare

# Make taxa uniquely identifiable 
tax_fix(
  ps.rare,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.rare

## Calculate alpha diversity ====
# Agglomerate to genus level
ps.rare.genus <- tax_glom(ps.rare, taxrank = "Genus", NArm=FALSE)

# Calculate
estimate_richness(ps.rare.genus, measures = c("Observed", "InvSimpson", "Shannon")) -> alphaDiv
### Add Samples column
alphaDiv$Samples <- rownames(alphaDiv)

# Calculate PD
pd(ps.rare.genus@otu_table, ps.rare.genus@phy_tree, include.root = TRUE) -> pd
### Add Samples column
pd$Samples <- rownames(pd)

### Combine pd with other metrics
merge(alphaDiv, pd, by = "Samples") -> alphaDiv

### Combine alpha metrics with sam_data
merge(as.matrix(ps.rare.genus@sam_data), alphaDiv, by = "Samples") -> alphaDiv

### Save dataframe
write.csv(alphaDiv, here::here("Data/07 - Alpha Diversity - Output/alphaDiv.csv"))

## Linear mixed effects models ====

# Specify "pretreatment"
alphaDiv$Pretreatment[alphaDiv$Time == "0"] <- "Pretreatment"
alphaDiv$PretreatmentDose[alphaDiv$Time == "0"] <- "Pretreatment"
alphaDiv$PretreatmentAbx[alphaDiv$Time == "0"] <- "Pretreatment"

# Set time as factor
alphaDiv$Time <- as.factor(alphaDiv$Time)

# Model: Affect of time, treatment, and the interaction on Shannon diversity
## Shannon
mod_shan <- lmer(Shannon ~ Treatment*Time + (1|Tank/CoralID), data = alphaDiv)
anova(mod_shan) -> anova_shan # Time and interaction of treatment and time is sig (p = 0.001991, p = 0.001571)
broom::tidy(anova_shan) -> anova_shan
anova_shan$Metric <- "Shannon"

# Check residuals of model
ggqqplot(residuals(mod_shan))

mod_shan_pw <- emmeans(mod_shan, pairwise ~ Time | Treatment, adjust = "tukey")
mod_shan_pw <- broom::tidy(mod_shan_pw$contrasts)
mod_shan_pw <- rstatix::add_significance(data = as.data.frame(mod_shan_pw), p.col = "adj.p.value")
mod_shan_pw$Metric <- "Shannon"

## Observed
mod_obs <- lmer(log(Observed) ~ Treatment*Time + (1|Tank/CoralID), data = alphaDiv)
anova(mod_obs) -> anova_obs # Time and interaction of treatment and time is sig (p = 0.035, p = 8.24e-05)
broom::tidy(anova_obs) -> anova_obs
anova_obs$Metric <- "Observed"
# Check residuals of model
ggqqplot(residuals(mod_obs))

mod_obs_pw <- emmeans(mod_obs, pairwise ~ Time | Treatment, adjust = "tukey")
mod_obs_pw <- broom::tidy(mod_obs_pw$contrasts)
mod_obs_pw <- rstatix::add_significance(data = as.data.frame(mod_obs_pw), p.col = "adj.p.value")
mod_obs_pw$Metric <- "Observed"

## Inverse Simpson
mod_invSimp <- lmer(log(InvSimpson) ~ Treatment*Time + (1|Tank/CoralID), data = alphaDiv)
anova(mod_invSimp) -> anova_invSimp # Time and interaction of treatment and time is sig (p = 1.48e-05, p = 0.0004)
broom::tidy(anova_invSimp) -> anova_invSimp
anova_invSimp$Metric <- "Inverse Simpson"
# Check residuals of model
ggqqplot(residuals(mod_invSimp))

mod_invSimp_pw <- emmeans(mod_invSimp, pairwise ~ Time | Treatment, adjust = "tukey")
mod_invSimp_pw <- broom::tidy(mod_invSimp_pw$contrasts)
mod_invSimp_pw <- rstatix::add_significance(data = as.data.frame(mod_invSimp_pw), p.col = "adj.p.value")
mod_invSimp_pw$Metric <- "Inverse Simpson"

## PD
mod_pd <- lmer(log(PD) ~ Treatment*Time + (1|Tank/CoralID), data = alphaDiv)
anova(mod_pd) -> anova_pd # Interaction only significant
broom::tidy(anova_pd) -> anova_pd
anova_pd$Metric <- "PD"
# Check residuals of model
ggqqplot(residuals(mod_pd))

mod_pd_pw <- emmeans(mod_pd, pairwise ~ Time | Treatment, adjust = "tukey")
mod_pd_pw <- broom::tidy(mod_pd_pw$contrasts)
mod_pd_pw <- rstatix::add_significance(data = as.data.frame(mod_pd_pw), p.col = "adj.p.value")
mod_pd_pw$Metric <- "PD"

# Combine stats output
## ANOVA output 
bind_rows(anova_shan, anova_obs, anova_invSimp, anova_pd) -> anova_output
write.csv(anova_output, here::here("Data/07 - Alpha Diversity - Output/anova_output.csv"))

## Pairwise output 
bind_rows(mod_shan_pw, mod_obs_pw, mod_invSimp_pw, mod_pd_pw) -> pairwise_stats

# Fix pairwise stats outputs to work with stat_pvalue_manual
pairwise_stats$group1 <- paste0(pairwise_stats$contrast)
pairwise_stats$group2 <- paste0(pairwise_stats$contrast)

pairwise_stats$group1[pairwise_stats$group1 == "Time0 - Time12"] <- "0"
pairwise_stats$group2[pairwise_stats$group2 == "Time0 - Time12"] <- "12"

pairwise_stats$group1[pairwise_stats$group1 == "Time0 - Time24"] <- "0"
pairwise_stats$group2[pairwise_stats$group2 == "Time0 - Time24"] <- "24"

pairwise_stats$group1[pairwise_stats$group1 == "Time0 - Time48"] <- "0"
pairwise_stats$group2[pairwise_stats$group2 == "Time0 - Time48"] <- "48"

pairwise_stats$group1[pairwise_stats$group1 == "Time0 - Time96"] <- "0"
pairwise_stats$group2[pairwise_stats$group2 == "Time0 - Time96"] <- "96"

pairwise_stats$group1[pairwise_stats$group1 == "Time12 - Time24"] <- "12"
pairwise_stats$group2[pairwise_stats$group2 == "Time12 - Time24"] <- "24"

pairwise_stats$group1[pairwise_stats$group1 == "Time12 - Time48"] <- "12"
pairwise_stats$group2[pairwise_stats$group2 == "Time12 - Time48"] <- "48"

pairwise_stats$group1[pairwise_stats$group1 == "Time12 - Time96"] <- "12"
pairwise_stats$group2[pairwise_stats$group2 == "Time12 - Time96"] <- "96"

pairwise_stats$group1[pairwise_stats$group1 == "Time24 - Time48"] <- "24"
pairwise_stats$group2[pairwise_stats$group2 == "Time24 - Time48"] <- "48"

pairwise_stats$group1[pairwise_stats$group1 == "Time24 - Time96"] <- "24"
pairwise_stats$group2[pairwise_stats$group2 == "Time24 - Time96"] <- "96"

pairwise_stats$group1[pairwise_stats$group1 == "Time48 - Time96"] <- "48"
pairwise_stats$group2[pairwise_stats$group2 == "Time48 - Time96"] <- "96"

## Add y.position column 
pairwise_stats$y.position <- 3.5

## Change column name
colnames(pairwise_stats)[10] <- "p.adj.signif"

## Plot Shannon
plot_shan <- ggplot(alphaDiv, aes(x = as.factor(Time), y = Shannon, color = Treatment)) + 
  scale_color_manual(values = c("#666666", "#FF7856", "#FF3D3D", "#9999FF","#3FA0FF", "#00BB00", "#BC71BB", "#006DDB", "#008600")) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA") + stat_summary(fun = mean, geom = "line", mapping = aes(group = Treatment, color = Treatment),
                                                               linewidth = 1.25, position = position_dodge(width = 0.9))
plot_shan <- plot_shan + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment, color = Treatment), size = 3,
                                      position = position_dodge(width = 0.9))
plot_shan <- plot_shan + geom_point(aes(color = Treatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_shan <- plot_shan + xlab("Time (Hours)") + ylab("Shannon") 
plot_shan <- plot_shan + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") + 
  stat_pvalue_manual(subset(pairwise_stats, Metric == "Shannon"), label = "p.adj.signif", hide.ns = TRUE) + 
  facet_wrap(~factor(Treatment, levels = c("Blank","Mixture Low", "Mixture High", "Ampicillin Low", "Streptomycin Low", 
                                           "Ciprofloxacin Low","Ampicillin High","Streptomycin High","Ciprofloxacin High"))) + theme(strip.text = element_text(face = "bold", size = 12))

# Change y.position for significance bars 
## Ampicillin low 
pairwise_stats$y.position[31] <- 3.2
pairwise_stats$y.position[32] <- 3.4
pairwise_stats$y.position[34] <- 3.6

## Streptomycin low
pairwise_stats$y.position[41] <- 3.2
pairwise_stats$y.position[42] <- 3.4
pairwise_stats$y.position[43] <- 3.6
# Rerun plot

ggplot2::ggsave(here::here("Data/07 - Alpha Diversity - Output/plot_shannon.png"), plot_shan,
                height = 450, width = 550, units = "mm",
                scale = 0.5, dpi = 1000)

## Plot observed
alphaDiv$Treatment <- factor(alphaDiv$Treatment, c("Blank", "Mixture Low", "Mixture High",
                                                   "Ampicillin Low", "Streptomycin Low", "Ciprofloxacin Low",
                                                   "Ampicillin High", "Streptomycin High", "Ciprofloxacin High"), 
                            order = TRUE)

plot_obs <- ggplot(alphaDiv, aes(x = as.factor(Time), y = Observed, color = Treatment)) + 
  scale_color_manual(values = c("#666666", "#FF7856", "#FF3D3D", "#9999FF","#3FA0FF", "#00BB00", "#BC71BB", "#006DDB", "#008600")) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA") + stat_summary(fun = mean, geom = "line", mapping = aes(group = Treatment, color = Treatment),
                                                               linewidth = 1.25, position = position_dodge(width = 0.9))
plot_obs <- plot_obs + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment, color = Treatment), size = 3,
                                      position = position_dodge(width = 0.9))
plot_obs <- plot_obs + geom_point(aes(color = Treatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_obs <- plot_obs + xlab("Time (Hours)") + ylab("Observed") 
plot_obs <- plot_obs + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") + 
  stat_pvalue_manual(subset(pairwise_stats, Metric == "Observed"), label = "p.adj.signif", hide.ns = TRUE) + 
  facet_wrap(~factor(Treatment, levels = c("Blank","Mixture Low", "Mixture High", "Ampicillin Low", "Streptomycin Low", "Ciprofloxacin Low","Ampicillin High","Streptomycin High","Ciprofloxacin High")))

# Change y.position for significance bars 
## Mixture High
pairwise_stats$y.position[144] <- 55
pairwise_stats$y.position[149] <- 50
## Ampicillin Low
pairwise_stats$y.position[103] <- 50
## Streptomycin Low
pairwise_stats$y.position[177] <- 55
pairwise_stats$y.position[179] <- 60
pairwise_stats$y.position[180] <- 65
## Ciprofloxacin High
pairwise_stats$y.position[123] <- 60
pairwise_stats$y.position[124] <- 65
# Rerun plot

## Plot inverse simpson
plot_invSimp <- ggplot(alphaDiv, aes(x = as.factor(Time), y = InvSimpson, color = Treatment)) + 
  scale_color_manual(values = c("#666666", "#FF7856", "#FF3D3D", "#9999FF","#3FA0FF", "#00BB00", "#BC71BB", "#006DDB", "#008600")) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA") + stat_summary(fun = mean, geom = "line", mapping = aes(group = Treatment, color = Treatment),
                                                               linewidth = 1.25, position = position_dodge(width = 0.9))
plot_invSimp <- plot_invSimp + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment, color = Treatment), size = 3,
                                    position = position_dodge(width = 0.9))
plot_invSimp <- plot_invSimp + geom_point(aes(color = Treatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_invSimp <- plot_invSimp + xlab("Time (Hours)") + ylab("Inverse Simpson") 
plot_invSimp <- plot_invSimp + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") + 
  stat_pvalue_manual(subset(pairwise_stats, Metric == "Inverse Simpson"), label = "p.adj.signif", hide.ns = TRUE) + 
  facet_wrap(~factor(Treatment, levels = c("Blank","Mixture Low", "Mixture High", "Ampicillin Low", "Streptomycin Low", "Ciprofloxacin Low","Ampicillin High","Streptomycin High","Ciprofloxacin High")))

# Change y.position for significance bars 
## Blank 
pairwise_stats$y.position[204] <- 25
## Mixture Low
pairwise_stats$y.position[242] <- 25
## Ampicillin Low
pairwise_stats$y.position[191] <- 20
pairwise_stats$y.position[192] <- 23
pairwise_stats$y.position[194] <- 25
## Streptomycin Low
pairwise_stats$y.position[261] <- 20
pairwise_stats$y.position[262] <- 23
pairwise_stats$y.position[263] <- 25
## Ciprofloxacin Low
pairwise_stats$y.position[228] <- 25
# Rerun plot

## Plot PD
plot_pd <- ggplot(alphaDiv, aes(x = as.factor(Time), y = PD, color = Treatment)) + 
  scale_color_manual(values = c("#666666", "#FF7856", "#FF3D3D", "#9999FF","#3FA0FF", "#00BB00", "#BC71BB", "#006DDB", "#008600")) +
  geom_boxplot(lwd = 1.1, outlier.color = "NA") + stat_summary(fun = mean, geom = "line", mapping = aes(group = Treatment, color = Treatment),
                                                               linewidth = 1.25, position = position_dodge(width = 0.9))
plot_pd <- plot_pd + stat_summary(fun = mean, geom = 'point', mapping = aes(group = Treatment, color = Treatment), size = 3,
                                            position = position_dodge(width = 0.9))
plot_pd <- plot_pd + geom_point(aes(color = Treatment), alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1)) + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))
plot_pd <- plot_pd + xlab("Time (Hours)") + ylab("Phylogenetic Diversity") 
plot_pd <- plot_pd + guides(color = guide_legend(title = "Treatment")) + theme(legend.position = "none") + 
  stat_pvalue_manual(subset(pairwise_stats, Metric == "PD"), label = "p.adj.signif", hide.ns = TRUE) + 
  facet_wrap(~factor(Treatment, levels = c("Blank","Mixture Low", "Mixture High", "Ampicillin Low", "Streptomycin Low", "Ciprofloxacin Low","Ampicillin High","Streptomycin High","Ciprofloxacin High")))

# Change y.position for significance bars 
## Mixture High
pairwise_stats$y.position[324] <- 6.5
pairwise_stats$y.position[330] <- 5.5
## Ampicillin Low
pairwise_stats$y.position[283] <- 5.5
## Ampicillin High
pairwise_stats$y.position[276] <- 5.5
# Rerun plot


# Save pairwise_stats output
write.csv(pairwise_stats, here::here("Data/07 - Alpha Diversity - Output/LMEM_pairwise.csv"))

# Save other alpha metric plots
ggarrange(plot_obs, plot_invSimp, plot_pd) -> alphaPlots

ggplot2::ggsave(here::here("Data/07 - Alpha Diversity - Output/plot_ObsSimpPD.png"), alphaPlots,
                height = 600, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)
