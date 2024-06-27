# Author: Sunni Patton
# Last edited: 06/03/24
# Title: Beta Diversity

## Set seed ====
set.seed(123)

## Load libraries ====
library(microViz)
library(pairwiseAdonis)
library(phyloseq)
library(ggpubr)
library(cowplot)
library(speedyseq)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(rstatix)

## Load data ====
ps.pruned <- readRDS(here::here("Data/05 - Phyloseq - Output/ps.pruned.RDS"))

## Edit sam_data ====
ps.pruned@sam_data$Pretreatment[ps.pruned@sam_data$Time == "0"] <- "Pretreatment"
ps.pruned@sam_data$PretreatmentAbx[ps.pruned@sam_data$Time == "0"] <- "Pretreatment"
ps.pruned@sam_data$PretreatmentDose[ps.pruned@sam_data$Time == "0"] <- "Pretreatment"

## Add additional columns (PretreamentTime, PretreatmentAbxTime, PretreatmentAbxDose)
ps.pruned@sam_data$PretreatmentTime <- paste0(ps.pruned@sam_data$Pretreatment, " ", ps.pruned@sam_data$Time)
ps.pruned@sam_data$PretreatmentTime[ps.pruned@sam_data$PretreatmentTime == "Pretreatment 0"] <- "Pretreatment"

ps.pruned@sam_data$PretreatmentAbxTime <- paste0(ps.pruned@sam_data$PretreatmentAbx, " ", ps.pruned@sam_data$Time)
ps.pruned@sam_data$PretreatmentAbxTime[ps.pruned@sam_data$PretreatmentAbxTime == "Pretreatment 0"] <- "Pretreatment"

ps.pruned@sam_data$PretreatmentDoseTime <- paste0(ps.pruned@sam_data$PretreatmentDose, " ", ps.pruned@sam_data$Time)
ps.pruned@sam_data$PretreatmentDoseTime[ps.pruned@sam_data$PretreatmentDoseTime == "Pretreatment 0"] <- "Pretreatment"

ps.pruned@sam_data$AntibioticTime <- paste0(ps.pruned@sam_data$Antibiotic, " ", ps.pruned@sam_data$Time)

## Validate phyloseq object and fix if necessary
ps.pruned <- phyloseq_validate(ps.pruned)

tax_fix(
  ps.pruned,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.pruned

## Are differences in beta diversity seen at T0? ====
subset_samples(ps.pruned, Time == "0") -> ps.pruned.T0

# rCLR transform data
ps.T0.rclr <- ps.pruned.T0 %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

# Calculate Euclidean distance
clr_euc_T0 <- phyloseq::distance(ps.T0.rclr, method = "euclidean")

# Extract sample data
sampledf_T0 <- data.frame(sample_data(ps.T0.rclr))

# adonis2 
adonis2(clr_euc_T0 ~ Antibiotic, data = sampledf_T0, perm = 999) # Not significant 
adonis2(clr_euc_T0 ~ Treatment, data = sampledf_T0, perm = 999) # Not significant 
adonis2(clr_euc_T0 ~ Tank, data = sampledf_T0, perm = 999) # Not significant

adonisT0 <- adonis2(clr_euc_T0 ~ Treatment, data = sampledf_T0, perm = 999) 
broom::tidy(adonisT0) -> adonisT0
adonisT0$Time <- "0"
adonisT0$Treatment <- "NA"

## Group T0 together as "Pretreatment" in PCA

## Are differences in beta diversity seen over time in blanks? ====
subset_samples(ps.pruned, Antibiotic == "Blank") -> ps.pruned.blank

# rCLR transform data
ps.blank.rclr <- ps.pruned.blank %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

# Calculate Euclidean distance
clr_euc_blank <- phyloseq::distance(ps.blank.rclr, method = "euclidean")

# Extract sample data
sampledf_blank <- data.frame(sample_data(ps.blank.rclr))

# adonis2 (does time affect non-treatment samples) 
adonis2(clr_euc_blank ~ Time, data = sampledf_blank, perm = 999) # Time is significant 

# pairwise adonis
pairwise.adonis(clr_euc_blank, sampledf_blank$Time, p.adjust.m = "fdr", perm = 999)

## Compare between groups within single time point ====
# Mixture 
subset_samples(ps.pruned, Antibiotic == "Blank" | Antibiotic == "Mixture") -> ps.blankmix

## rCLR transform data
ps.blankmix.rclr <- ps.blankmix %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate euclidean distance 
### T12 
clr_euc_blankmix_T12 <- phyloseq::distance(subset_samples(ps.blankmix.rclr, Time == "12"), method = "euclidean")
sampledf_blankmix_T12 <- data.frame(sample_data(subset_samples(ps.blankmix.rclr, Time == "12")))

### T24 
clr_euc_blankmix_T24 <- phyloseq::distance(subset_samples(ps.blankmix.rclr, Time == "24"), method = "euclidean")
sampledf_blankmix_T24 <- data.frame(sample_data(subset_samples(ps.blankmix.rclr, Time == "24")))

### T48 
clr_euc_blankmix_T48 <- phyloseq::distance(subset_samples(ps.blankmix.rclr, Time == "48"), method = "euclidean")
sampledf_blankmix_T48 <- data.frame(sample_data(subset_samples(ps.blankmix.rclr, Time == "48")))

### T96 
clr_euc_blankmix_T96 <- phyloseq::distance(subset_samples(ps.blankmix.rclr, Time == "96"), method = "euclidean")
sampledf_blankmix_T96 <- data.frame(sample_data(subset_samples(ps.blankmix.rclr, Time == "96")))

## adonis2
### T12
adonis_blankmix_T12 <- adonis2(clr_euc_blankmix_T12 ~ Treatment, data = sampledf_blankmix_T12, perm = 999) # Sig
adonis_blankmix_T12 <- broom::tidy(adonis_blankmix_T12)
adonis_blankmix_T12$Time <- "12"
adonis_blankmix_T12$Treatment <- "Blank-Mixture High-Mixture Low"

### T24
adonis_blankmix_T24 <- adonis2(clr_euc_blankmix_T24 ~ Treatment, data = sampledf_blankmix_T24, perm = 999) # Treatment significant 
adonis_blankmix_T24 <- broom::tidy(adonis_blankmix_T24)
adonis_blankmix_T24$Time <- "24"
adonis_blankmix_T24$Treatment <- "Blank-Mixture High-Mixture Low"

### T48
adonis_blankmix_T48 <- adonis2(clr_euc_blankmix_T48 ~ Treatment, data = sampledf_blankmix_T48, perm = 999) # Treatment significant
adonis_blankmix_T48 <- broom::tidy(adonis_blankmix_T48)
adonis_blankmix_T48$Time <- "48"
adonis_blankmix_T48$Treatment <- "Blank-Mixture High-Mixture Low"

### T96
adonis_blankmix_T96 <- adonis2(clr_euc_blankmix_T96 ~ Treatment, data = sampledf_blankmix_T96, perm = 999) # Treatment significant
adonis_blankmix_T96 <- broom::tidy(adonis_blankmix_T96)
adonis_blankmix_T96$Time <- "96"
adonis_blankmix_T96$Treatment <- "Blank-Mixture High-Mixture Low"

## Combine adonis2 output
adonis_blankmix <- rbind(adonis_blankmix_T12, adonis_blankmix_T24, adonis_blankmix_T48, adonis_blankmix_T96)

## Pairwise adonis
pwadonis_blankmix_T12 <- pairwise.adonis(clr_euc_blankmix_T12, sampledf_blankmix_T12$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between mix high and blank, mix low and blank, but not sig between mix high and low
pwadonis_blankmix_T12$Time <- "12"

pwadonis_blankmix_T24 <- pairwise.adonis(clr_euc_blankmix_T24, sampledf_blankmix_T24$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between mix high and blank, mix low and blank, but not sig between mix high and low
pwadonis_blankmix_T24$Time <- "24"

pwadonis_blankmix_T48 <- pairwise.adonis(clr_euc_blankmix_T48, sampledf_blankmix_T48$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between mix high and blank, mix low and blank, but not sig between mix high and low
pwadonis_blankmix_T48$Time <- "48"

pwadonis_blankmix_T96 <- pairwise.adonis(clr_euc_blankmix_T96, sampledf_blankmix_T96$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between mix high and blank, mix low and blank, and between mix high and mix low
pwadonis_blankmix_T96$Time <- "96"

## Combine pairwise adonis dataframes 
pwadonis_blankmix <- rbind(pwadonis_blankmix_T12, pwadonis_blankmix_T24, pwadonis_blankmix_T48, pwadonis_blankmix_T96)

# Ampicillin
subset_samples(ps.pruned, Antibiotic == "Blank" | Antibiotic == "Ampicillin") -> ps.blankamp

## rCLR transform data
ps.blankamp.rclr <- ps.blankamp %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate euclidean distance 
### T12 Antibiotic
clr_euc_blankamp_T12 <- phyloseq::distance(subset_samples(ps.blankamp.rclr, Time == "12"), method = "euclidean")
sampledf_blankamp_T12 <- data.frame(sample_data(subset_samples(ps.blankamp.rclr, Time == "12")))

### T24 
clr_euc_blankamp_T24 <- phyloseq::distance(subset_samples(ps.blankamp.rclr, Time == "24"), method = "euclidean")
sampledf_blankamp_T24 <- data.frame(sample_data(subset_samples(ps.blankamp.rclr, Time == "24")))

### T48 
clr_euc_blankamp_T48 <- phyloseq::distance(subset_samples(ps.blankamp.rclr, Time == "48"), method = "euclidean")
sampledf_blankamp_T48 <- data.frame(sample_data(subset_samples(ps.blankamp.rclr, Time == "48")))

### T96 
clr_euc_blankamp_T96 <- phyloseq::distance(subset_samples(ps.blankamp.rclr, Time == "96"), method = "euclidean")
sampledf_blankamp_T96 <- data.frame(sample_data(subset_samples(ps.blankamp.rclr, Time == "96")))

## adonis2
### T12 
adonis_blankamp_T12 <- adonis2(clr_euc_blankamp_T12 ~ Treatment, data = sampledf_blankamp_T12, perm = 999) # Treatment significant
adonis_blankamp_T12 <- broom::tidy(adonis_blankamp_T12)
adonis_blankamp_T12$Time <- "12"
adonis_blankamp_T12$Treatment <- "Blank-Ampicillin High-Ampicillin Low"

### T24 Treatment
adonis_blankamp_T24 <- adonis2(clr_euc_blankamp_T24 ~ Treatment, data = sampledf_blankamp_T24, perm = 999) # Treatment significant 
adonis_blankamp_T24 <- broom::tidy(adonis_blankamp_T24)
adonis_blankamp_T24$Time <- "24"
adonis_blankamp_T24$Treatment <- "Blank-Ampicillin High-Ampicillin Low"

### T48 Treatment
adonis_blankamp_T48 <- adonis2(clr_euc_blankamp_T48 ~ Treatment, data = sampledf_blankamp_T48, perm = 999) # Treatment significant
adonis_blankamp_T48 <- broom::tidy(adonis_blankamp_T48)
adonis_blankamp_T48$Time <- "48"
adonis_blankamp_T48$Treatment <- "Blank-Ampicillin High-Ampicillin Low"

### T96 Treatment
adonis_blankamp_T96 <- adonis2(clr_euc_blankamp_T96 ~ Treatment, data = sampledf_blankamp_T96, perm = 999) # Treatment significant
adonis_blankamp_T96 <- broom::tidy(adonis_blankamp_T96)
adonis_blankamp_T96$Time <- "96"
adonis_blankamp_T96$Treatment <- "Blank-Ampicillin High-Ampicillin Low"

### Combine adonis2 output
adonis_blankamp <- rbind(adonis_blankamp_T12, adonis_blankamp_T24, adonis_blankamp_T48, adonis_blankamp_T96)

# pairwise adonis
pwadonis_blankamp_T12 <- pairwise.adonis(clr_euc_blankamp_T12, sampledf_blankamp_T12$Treatment, p.adjust.m = "fdr", perm = 999) # No longer significant after p value correction
pwadonis_blankamp_T12$Time <- "12"

pwadonis_blankamp_T24 <- pairwise.adonis(clr_euc_blankamp_T24, sampledf_blankamp_T24$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between amp high and blank, amp low and blank, but not sig between amp high and low
pwadonis_blankamp_T24$Time <- "24"

pwadonis_blankamp_T48 <- pairwise.adonis(clr_euc_blankamp_T48, sampledf_blankamp_T48$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between amp high and blank, amp low and blank, and between amp high and amp low
pwadonis_blankamp_T48$Time <- "48"

pwadonis_blankamp_T96 <- pairwise.adonis(clr_euc_blankamp_T96, sampledf_blankamp_T96$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between amp high and blank, amp low and blank, and between amp high and amp low
pwadonis_blankamp_T96$Time <- "96"

## Combine pairwise adonis dataframes
pwadonis_blankamp <- rbind(pwadonis_blankamp_T12, pwadonis_blankamp_T24, pwadonis_blankamp_T48, pwadonis_blankamp_T96)

## Streptomycin
subset_samples(ps.pruned, Antibiotic == "Blank" | Antibiotic == "Streptomycin") -> ps.blankstrep

## rCLR transform data
ps.blankstrep.rclr <- ps.blankstrep %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate euclidean distance 
### T12 
clr_euc_blankstrep_T12 <- phyloseq::distance(subset_samples(ps.blankstrep.rclr, Time == "12"), method = "euclidean")
sampledf_blankstrep_T12 <- data.frame(sample_data(subset_samples(ps.blankstrep.rclr, Time == "12")))

### T24 
clr_euc_blankstrep_T24 <- phyloseq::distance(subset_samples(ps.blankstrep.rclr, Time == "24"), method = "euclidean")
sampledf_blankstrep_T24 <- data.frame(sample_data(subset_samples(ps.blankstrep.rclr, Time == "24")))

### T48 
clr_euc_blankstrep_T48 <- phyloseq::distance(subset_samples(ps.blankstrep.rclr, Time == "48"), method = "euclidean")
sampledf_blankstrep_T48 <- data.frame(sample_data(subset_samples(ps.blankstrep.rclr, Time == "48")))

### T96 
clr_euc_blankstrep_T96 <- phyloseq::distance(subset_samples(ps.blankstrep.rclr, Time == "96"), method = "euclidean")
sampledf_blankstrep_T96 <- data.frame(sample_data(subset_samples(ps.blankstrep.rclr, Time == "96")))

## adonis2
### T12 
adonis_blankstrep_T12 <- adonis2(clr_euc_blankstrep_T12 ~ Treatment, data = sampledf_blankstrep_T12, perm = 999) # Treatment significant
adonis_blankstrep_T12 <- broom::tidy(adonis_blankstrep_T12)
adonis_blankstrep_T12$Time <- "12"
adonis_blankstrep_T12$Treatment <- "Blank-Streptomycin High-Streptomycin Low"

### T24 
adonis_blankstrep_T24 <- adonis2(clr_euc_blankstrep_T24 ~ Treatment, data = sampledf_blankstrep_T24, perm = 999) # Treatment significant 
adonis_blankstrep_T24 <- broom::tidy(adonis_blankstrep_T24)
adonis_blankstrep_T24$Time <- "24"
adonis_blankstrep_T24$Treatment <- "Blank-Streptomycin High-Streptomycin Low"

### T48 
adonis_blankstrep_T48 <- adonis2(clr_euc_blankstrep_T48 ~ Treatment, data = sampledf_blankstrep_T48, perm = 999) # Treatment significant
adonis_blankstrep_T48 <- broom::tidy(adonis_blankstrep_T48)
adonis_blankstrep_T48$Time <- "48"
adonis_blankstrep_T48$Treatment <- "Blank-Streptomycin High-Streptomycin Low"

### T96 
adonis_blankstrep_T96 <- adonis2(clr_euc_blankstrep_T96 ~ Treatment, data = sampledf_blankstrep_T96, perm = 999) # Treatment significant
adonis_blankstrep_T96 <- broom::tidy(adonis_blankstrep_T96)
adonis_blankstrep_T96$Time <- "96"
adonis_blankstrep_T96$Treatment <- "Blank-Streptomycin High-Streptomycin Low"

## Combine adonis2 dataframes
adonis_blankstrep <- rbind(adonis_blankstrep_T12, adonis_blankstrep_T24, adonis_blankstrep_T48, adonis_blankstrep_T96)

## Pairwise adonis
pwadonis_blankstrep_T12 <- pairwise.adonis(clr_euc_blankstrep_T12, sampledf_blankstrep_T12$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between strep high and blank, strep low and blank, but not sig between strep high and low
pwadonis_blankstrep_T12$Time <- "12"

pwadonis_blankstrep_T24 <- pairwise.adonis(clr_euc_blankstrep_T24, sampledf_blankstrep_T24$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between strep high and blank, strep low and blank, and between strep high and strep low
pwadonis_blankstrep_T24$Time <- "24"

pwadonis_blankstrep_T48 <- pairwise.adonis(clr_euc_blankstrep_T48, sampledf_blankstrep_T48$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between strep high and blank, strep low and blank, and between strep high and strep low
pwadonis_blankstrep_T48$Time <- "48"

pwadonis_blankstrep_T96 <- pairwise.adonis(clr_euc_blankstrep_T96, sampledf_blankstrep_T96$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between strep high and blank, strep low and blank, and between strep high and strep low
pwadonis_blankstrep_T96$Time <- "96"

## Combine pairwise adonis dataframes
pwadonis_blankstrep <- rbind(pwadonis_blankstrep_T12, pwadonis_blankstrep_T24, pwadonis_blankstrep_T48, pwadonis_blankstrep_T96)

## Ciprofloxacin
subset_samples(ps.pruned, Antibiotic == "Blank" | Antibiotic == "Ciprofloxacin") -> ps.blankcip

## rCLR transform data
ps.blankcip.rclr <- ps.blankcip %>%
  tax_transform(trans = "rclr", rank = "Genus") %>%
  ord_calc(method = "PCA")

## Calculate euclidean distance 
### T12 
clr_euc_blankcip_T12 <- phyloseq::distance(subset_samples(ps.blankcip.rclr, Time == "12"), method = "euclidean")
sampledf_blankcip_T12 <- data.frame(sample_data(subset_samples(ps.blankcip.rclr, Time == "12")))

### T24 
clr_euc_blankcip_T24 <- phyloseq::distance(subset_samples(ps.blankcip.rclr, Time == "24"), method = "euclidean")
sampledf_blankcip_T24 <- data.frame(sample_data(subset_samples(ps.blankcip.rclr, Time == "24")))

### T48 
clr_euc_blankcip_T48 <- phyloseq::distance(subset_samples(ps.blankcip.rclr, Time == "48"), method = "euclidean")
sampledf_blankcip_T48 <- data.frame(sample_data(subset_samples(ps.blankcip.rclr, Time == "48")))

### T96 
clr_euc_blankcip_T96 <- phyloseq::distance(subset_samples(ps.blankcip.rclr, Time == "96"), method = "euclidean")
sampledf_blankcip_T96 <- data.frame(sample_data(subset_samples(ps.blankcip.rclr, Time == "96")))

## adonis2
### T12 
adonis_blankcip_T12 <- adonis2(clr_euc_blankcip_T12 ~ Treatment, data = sampledf_blankcip_T12, perm = 999) # Treatment significant
adonis_blankcip_T12 <- broom::tidy(adonis_blankcip_T12)
adonis_blankcip_T12$Time <- "12"
adonis_blankcip_T12$Treatment <- "Blank-Ciprofloxacin High-Ciprofloxacin Low"

### T24 
adonis_blankcip_T24 <- adonis2(clr_euc_blankcip_T24 ~ Treatment, data = sampledf_blankcip_T24, perm = 999) # Treatment significant 
adonis_blankcip_T24 <- broom::tidy(adonis_blankcip_T24)
adonis_blankcip_T24$Time <- "24"
adonis_blankcip_T24$Treatment <- "Blank-Ciprofloxacin High-Ciprofloxacin Low"

### T48 
adonis_blankcip_T48 <- adonis2(clr_euc_blankcip_T48 ~ Treatment, data = sampledf_blankcip_T48, perm = 999) # Treatment significant
adonis_blankcip_T48 <- broom::tidy(adonis_blankcip_T48)
adonis_blankcip_T48$Time <- "48"
adonis_blankcip_T48$Treatment <- "Blank-Ciprofloxacin High-Ciprofloxacin Low"

### T96 
adonis_blankcip_T96 <- adonis2(clr_euc_blankcip_T96 ~ Treatment, data = sampledf_blankcip_T96, perm = 999) # Treatment significant
adonis_blankcip_T96 <- broom::tidy(adonis_blankcip_T96)
adonis_blankcip_T96$Time <- "96"
adonis_blankcip_T96$Treatment <- "Blank-Ciprofloxacin High-Ciprofloxacin Low"

## Combine adonis2 output
adonis_blankcip <- rbind(adonis_blankcip_T12, adonis_blankcip_T24, adonis_blankcip_T48, adonis_blankcip_T96)

# pairwise adonis
pwadonis_blankcip_T12 <- pairwise.adonis(clr_euc_blankcip_T12, sampledf_blankcip_T12$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between strep high and blank, strep low and blank, but not sig between strep high and low
pwadonis_blankcip_T12$Time <- "12"

pwadonis_blankcip_T24 <- pairwise.adonis(clr_euc_blankcip_T24, sampledf_blankcip_T24$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between strep high and blank, strep low and blank, and between strep high and strep low
pwadonis_blankcip_T24$Time <- "24"

pwadonis_blankcip_T48 <- pairwise.adonis(clr_euc_blankcip_T48, sampledf_blankcip_T48$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between strep high and blank, strep low and blank, and between strep high and strep low
pwadonis_blankcip_T48$Time <- "48"

pwadonis_blankcip_T96 <- pairwise.adonis(clr_euc_blankcip_T96, sampledf_blankcip_T96$Treatment, p.adjust.m = "fdr", perm = 999) # Sig between strep high and blank, strep low and blank, and between strep high and strep low
pwadonis_blankcip_T96$Time <- "96"

## Combine pairwise adonis
pwadonis_blankcip <- rbind(pwadonis_blankcip_T12, pwadonis_blankcip_T24, pwadonis_blankcip_T48, pwadonis_blankcip_T96)

## Combine all adonis and pairwise adonis outputs ====
adonis_output <- rbind(adonisT0, adonis_blankmix, adonis_blankamp, adonis_blankstrep, adonis_blankcip)

write.csv(adonis_output, here::here("Data/08 - Beta Diversity - Output/02 - adonis_output_treat.csv"))

pwadonis_output <- rbind(pwadonis_blankmix, pwadonis_blankamp, pwadonis_blankstrep, pwadonis_blankcip)

write.csv(pwadonis_output, here::here("Data/08 - Beta Diversity - Output/03 - pwadonis_output.csv"))

## Plot Mixture + Blank ====
# Set order
ps.blankmix.rclr@sam_data$Pretreatment <- factor(ps.blankmix.rclr@sam_data$Pretreatment, c("Pretreatment", "Blank", "Mixture Low", "Mixture High"), 
                                                 order = TRUE)
PCA_blankmix <- ps.blankmix.rclr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + coord_fixed(ratio = 1, clip = "off", xlim = c(-2, 2)) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1) + facet_wrap(~Time, ncol = 5, nrow = 1) + theme_bw() + 
  scale_color_manual(values = c("#000000","#666666","#FF7856", "#FF3D3D")) + theme_bw() + 
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold")) + guides(color = guide_legend(title = "Treatment")) + 
  theme(strip.text = element_text(face = "bold", size = 12)) 


## Plot Ampicillin + Blank ====
# Set order
ps.blankamp.rclr@sam_data$Pretreatment <- factor(ps.blankamp.rclr@sam_data$Pretreatment, c("Pretreatment", "Blank", "Ampicillin Low", "Ampicillin High"), 
                                                 order = TRUE)

PCA_blankamp <- ps.blankamp.rclr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + coord_fixed(ratio = 1, clip = "off", xlim = c(-2, 2)) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1) + facet_wrap(~Time, ncol = 5, nrow = 1) + theme_bw() + 
  scale_color_manual(values = c("#000000","#666666","#9999FF", "#BC71BB")) + theme_bw() + 
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold")) + guides(color = guide_legend(title = "Treatment")) + 
  theme(strip.text = element_text(face = "bold", size = 12))

## Plot Streptomycin + Blank ====
# Set order
ps.blankstrep.rclr@sam_data$Pretreatment <- factor(ps.blankstrep.rclr@sam_data$Pretreatment, c("Pretreatment", "Blank", "Streptomycin Low", "Streptomycin High"), 
                                                 order = TRUE)

PCA_blankstrep <- ps.blankstrep.rclr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + coord_fixed(ratio = 1, clip = "off", xlim = c(-2, 2)) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1) + facet_wrap(~Time, ncol = 5, nrow = 1) + theme_bw() + 
  scale_color_manual(values = c("#000000","#666666","#3FA0FF", "#006DDB")) + theme_bw() + 
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold")) + guides(color = guide_legend(title = "Treatment")) + 
  theme(strip.text = element_text(face = "bold", size = 12)) 

## Plot Ciprofloxacin + Blank ====
# Set order
ps.blankcip.rclr@sam_data$Pretreatment <- factor(ps.blankcip.rclr@sam_data$Pretreatment, c("Pretreatment", "Blank", "Ciprofloxacin Low", "Ciprofloxacin High"), 
                                                   order = TRUE)

PCA_blankcip <- ps.blankcip.rclr %>%
  ord_plot( 
    colour = "Pretreatment",
    plot_taxa = FALSE, 
    auto_caption = NA
  ) + coord_fixed(ratio = 1, clip = "off", xlim = c(-2, 2)) + stat_ellipse(aes(group = Pretreatment, color = Pretreatment), linewidth = 1) + facet_wrap(~Time, ncol = 5, nrow = 1) + theme_bw() + 
  scale_color_manual(values = c("#000000","#666666","#00BB00", "#008600")) + theme_bw() + 
  theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), 
        title = element_text(face = "bold")) + guides(color = guide_legend(title = "Treatment")) + 
  theme(strip.text = element_text(face = "bold", size = 12))

## Beta dispersion T0 ====
disp_T0 <- betadisper(clr_euc_T0, sampledf_T0$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_T0) # Not sig p = 0.2125

stats_T0 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_T0, sampledf_T0$Treatment, type = "centroid", 
                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Significant between Blank-Mixture High

stats_T0$Time <- "0"

# Change column name from x to p.adj
colnames(stats_T0)[colnames(stats_T0) == 'x'] <- "p.adj"

# Add significance column
stats_T0 <- add_significance(stats_T0, p.col = "p.adj")

stats_T0$Group <- "Pretreatment"

df_disp_T0 <- data.frame(disp_T0$distances)
colnames(df_disp_T0)[colnames(df_disp_T0) == 'disp_T0.distances'] <- "Distances"
df_disp_T0$Time <- 0
df_disp_T0$Group <- "Pretreatment"


## Beta dispersion Mixture + Blank ====
## T12
disp_blankmix_T12 <- betadisper(clr_euc_blankmix_T12, sampledf_blankmix_T12$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankmix_T12) # Not sig p = 0.9791

stats_blankmix_T12 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankmix_T12, sampledf_blankmix_T12$Treatment, type = "centroid", 
                              bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_blankmix_T12$Time <- "12"


df_disp_blankmix_T12 <- data.frame(disp_blankmix_T12$distances)
colnames(df_disp_blankmix_T12)[colnames(df_disp_blankmix_T12) == 'disp_blankmix_T12.distances'] <- "Distances"
df_disp_blankmix_T12$Time <- 12

## T24
disp_blankmix_T24 <- betadisper(clr_euc_blankmix_T24, sampledf_blankmix_T24$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankmix_T24) # Not sig p = 0.5139

stats_blankmix_T24 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankmix_T24, sampledf_blankmix_T24$Treatment, type = "centroid", 
                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_blankmix_T24$Time <- "24"

df_disp_blankmix_T24 <- data.frame(disp_blankmix_T24$distances)
colnames(df_disp_blankmix_T24)[colnames(df_disp_blankmix_T24) == 'disp_blankmix_T24.distances'] <- "Distances"
df_disp_blankmix_T24$Time <- 24

## T48
disp_blankmix_T48 <- betadisper(clr_euc_blankmix_T48, sampledf_blankmix_T48$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankmix_T48) # Sig p = 0.009881

stats_blankmix_T48 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankmix_T48, sampledf_blankmix_T48$Treatment, type = "centroid", 
                              bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Significant between Blank-Mixture High

stats_blankmix_T48$Time <- "48"

df_disp_blankmix_T48 <- data.frame(disp_blankmix_T48$distances)
colnames(df_disp_blankmix_T48)[colnames(df_disp_blankmix_T48) == 'disp_blankmix_T48.distances'] <- "Distances"
df_disp_blankmix_T48$Time <- 48

## T96
disp_blankmix_T96 <- betadisper(clr_euc_blankmix_T96, sampledf_blankmix_T96$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankmix_T96) # Not sig p = 0.2022

stats_blankmix_T96 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankmix_T96, sampledf_blankmix_T96$Treatment, type = "centroid", 
                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Significant between Blank-Mixture High

stats_blankmix_T96$Time <- "96"

df_disp_blankmix_T96 <- data.frame(disp_blankmix_T96$distances)
colnames(df_disp_blankmix_T96)[colnames(df_disp_blankmix_T96) == 'disp_blankmix_T96.distances'] <- "Distances"
df_disp_blankmix_T96$Time <- 96


## Combine Blank Mix distance dataframes 
df_disp_blankmix <- rbind(df_disp_blankmix_T12, df_disp_blankmix_T24, df_disp_blankmix_T48, df_disp_blankmix_T96)
df_disp_blankmix$Group <- "Mixture"

## Combine Blank Mix stats dataframes
stats_blankmix <- rbind(stats_blankmix_T12, stats_blankmix_T24, stats_blankmix_T48, stats_blankmix_T96)

# Change column name from x to p.adj
colnames(stats_blankmix)[colnames(stats_blankmix) == 'x'] <- "p.adj"

# Add significance column
stats_blankmix <- add_significance(stats_blankmix, p.col = "p.adj")

# Add group column
stats_blankmix$Group <- "Mixture"


## Beta dispersion Ampicillin + Blank ====
## T12
disp_blankamp_T12 <- betadisper(clr_euc_blankamp_T12, sampledf_blankamp_T12$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankamp_T12) # Not sig p = 0.1962

stats_blankamp_T12 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankamp_T12, sampledf_blankamp_T12$Treatment, type = "centroid", 
                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_blankamp_T12$Time <- "12"

df_disp_blankamp_T12 <- data.frame(disp_blankamp_T12$distances)
colnames(df_disp_blankamp_T12)[colnames(df_disp_blankamp_T12) == 'disp_blankamp_T12.distances'] <- "Distances"
df_disp_blankamp_T12$Time <- 12

## T24
disp_blankamp_T24 <- betadisper(clr_euc_blankamp_T24, sampledf_blankamp_T24$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankamp_T24) # Not sig p = 0.6054

stats_blankamp_T24 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankamp_T24, sampledf_blankamp_T24$Treatment, type = "centroid", 
                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))

stats_blankamp_T24$Time <- 24

df_disp_blankamp_T24 <- data.frame(disp_blankamp_T24$distances)
colnames(df_disp_blankamp_T24)[colnames(df_disp_blankamp_T24) == 'disp_blankamp_T24.distances'] <- "Distances"
df_disp_blankamp_T24$Time <- 24

## T48
disp_blankamp_T48 <- betadisper(clr_euc_blankamp_T48, sampledf_blankamp_T48$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankamp_T48) # Sig p = 0.01869 

stats_blankamp_T48 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankamp_T48, sampledf_blankamp_T48$Treatment, type = "centroid", 
                              bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Significant between amp high and amp low

stats_blankamp_T48$Time <- 48

df_disp_blankamp_T48 <- data.frame(disp_blankamp_T48$distances)
colnames(df_disp_blankamp_T48)[colnames(df_disp_blankamp_T48) == 'disp_blankamp_T48.distances'] <- "Distances"
df_disp_blankamp_T48$Time <- 48

## T96
disp_blankamp_T96 <- betadisper(clr_euc_blankamp_T96, sampledf_blankamp_T96$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankamp_T96) # Not sig p = 0.317

stats_blankamp_T96 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankamp_T96, sampledf_blankamp_T96$Treatment, type = "centroid", 
                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Significant between amp high and amp low

stats_blankamp_T96$Time <- 96


df_disp_blankamp_T96 <- data.frame(disp_blankamp_T96$distances)
colnames(df_disp_blankamp_T96)[colnames(df_disp_blankamp_T96) == 'disp_blankamp_T96.distances'] <- "Distances"
df_disp_blankamp_T96$Time <- 96

## Combine Blank Amp distance dataframes 
df_disp_blankamp <- rbind(df_disp_blankamp_T12, df_disp_blankamp_T24, df_disp_blankamp_T48, df_disp_blankamp_T96)
df_disp_blankamp$Group <- "Ampicillin"

## Combine Blank Amp stats dataframes
stats_blankamp <- rbind(stats_blankamp_T12, stats_blankamp_T24, stats_blankamp_T48, stats_blankamp_T96)

# Change column name from x to p.adj
colnames(stats_blankamp)[colnames(stats_blankamp) == 'x'] <- "p.adj"

# Add significance column
stats_blankamp <- add_significance(stats_blankamp, p.col = "p.adj")

# Add group column
stats_blankamp$Group <- "Ampicillin"

## Beta dispersion Streptomycin + Blank ====
## T12
disp_blankstrep_T12 <- betadisper(clr_euc_blankstrep_T12, sampledf_blankstrep_T12$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankstrep_T12) # Not sig p = 0.1268

stats_blankstrep_T12 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankstrep_T12, sampledf_blankstrep_T12$Treatment, type = "centroid", 
                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_blankstrep_T12$Time <- "12"

df_disp_blankstrep_T12 <- data.frame(disp_blankstrep_T12$distances)
colnames(df_disp_blankstrep_T12)[colnames(df_disp_blankstrep_T12) == 'disp_blankstrep_T12.distances'] <- "Distances"
df_disp_blankstrep_T12$Time <- 12

## T24
disp_blankstrep_T24 <- betadisper(clr_euc_blankstrep_T24, sampledf_blankstrep_T24$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankstrep_T24) # Not sig p = 0.9639

stats_blankstrep_T24 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankstrep_T24, sampledf_blankstrep_T24$Treatment, type = "centroid", 
                                                                  bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_blankstrep_T24$Time <- "24"

df_disp_blankstrep_T24 <- data.frame(disp_blankstrep_T24$distances)
colnames(df_disp_blankstrep_T24)[colnames(df_disp_blankstrep_T24) == 'disp_blankstrep_T24.distances'] <- "Distances"
df_disp_blankstrep_T24$Time <- 24

## T48
disp_blankstrep_T48 <- betadisper(clr_euc_blankstrep_T48, sampledf_blankstrep_T48$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankstrep_T48) # Not sig p = 0.7286 

stats_blankstrep_T48 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankstrep_T48, sampledf_blankstrep_T48$Treatment, type = "centroid", 
                                                                  bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_blankstrep_T48$Time <- "48"

df_disp_blankstrep_T48 <- data.frame(disp_blankstrep_T48$distances)
colnames(df_disp_blankstrep_T48)[colnames(df_disp_blankstrep_T48) == 'disp_blankstrep_T48.distances'] <- "Distances"
df_disp_blankstrep_T48$Time <- 48

## T96
disp_blankstrep_T96 <- betadisper(clr_euc_blankstrep_T96, sampledf_blankstrep_T96$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankstrep_T96) # Not sig p = 0.1501

stats_blankstrep_T96 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankstrep_T96, sampledf_blankstrep_T96$Treatment, type = "centroid", 
                                                                  bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr'))
stats_blankstrep_T96$Time <- "96"

df_disp_blankstrep_T96 <- data.frame(disp_blankstrep_T96$distances)
colnames(df_disp_blankstrep_T96)[colnames(df_disp_blankstrep_T96) == 'disp_blankstrep_T96.distances'] <- "Distances"
df_disp_blankstrep_T96$Time <- 96

## Combine Blank Strep distance dataframes 
df_disp_blankstrep <- rbind(df_disp_blankstrep_T12, df_disp_blankstrep_T24, df_disp_blankstrep_T48, df_disp_blankstrep_T96)
df_disp_blankstrep$Group <- "Streptomycin"

## Combine Blank Strep stats dataframes
stats_blankstrep <- rbind(stats_blankstrep_T12, stats_blankstrep_T24, stats_blankstrep_T48, stats_blankstrep_T96)

# Change column name from x to p.adj
colnames(stats_blankstrep)[colnames(stats_blankstrep) == 'x'] <- "p.adj"

# Add significance column
stats_blankstrep <- add_significance(stats_blankstrep, p.col = "p.adj")

# Add group column
stats_blankstrep$Group <- "Streptomycin"

## Beta dispersion Ciprofloxacin + Blank ====
## T12
disp_blankcip_T12 <- betadisper(clr_euc_blankcip_T12, sampledf_blankcip_T12$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankcip_T12) # Sig p = 0.04421 

stats_blankcip_T12 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankcip_T12, sampledf_blankcip_T12$Treatment, type = "centroid", 
                              bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # No longer significant after correction

stats_blankcip_T12$Time <- 12

df_disp_blankcip_T12 <- data.frame(disp_blankcip_T12$distances)
colnames(df_disp_blankcip_T12)[colnames(df_disp_blankcip_T12) == 'disp_blankcip_T12.distances'] <- "Distances"
df_disp_blankcip_T12$Time <- 12

## T24
disp_blankcip_T24 <- betadisper(clr_euc_blankcip_T24, sampledf_blankcip_T24$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankcip_T24) # Sig p = 0.0194 

stats_blankcip_T24 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankcip_T24, sampledf_blankcip_T24$Treatment, type = "centroid", 
                              bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Significant between blank-cip low and cip high-cip low

stats_blankcip_T24$Time <- 24

df_disp_blankcip_T24 <- data.frame(disp_blankcip_T24$distances)
colnames(df_disp_blankcip_T24)[colnames(df_disp_blankcip_T24) == 'disp_blankcip_T24.distances'] <- "Distances"
df_disp_blankcip_T24$Time <- 24

## T48
disp_blankcip_T48 <- betadisper(clr_euc_blankcip_T48, sampledf_blankcip_T48$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankcip_T48) # Not sig p = 0.05338  

stats_blankcip_T48 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankcip_T48, sampledf_blankcip_T48$Treatment, type = "centroid", 
                                                                bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Significant between blank-cip low and cip high-cip low

stats_blankcip_T48$Time <- 48

df_disp_blankcip_T48 <- data.frame(disp_blankcip_T48$distances)
colnames(df_disp_blankcip_T48)[colnames(df_disp_blankcip_T48) == 'disp_blankcip_T48.distances'] <- "Distances"
df_disp_blankcip_T48$Time <- 48

## T96
disp_blankcip_T96 <- betadisper(clr_euc_blankcip_T96, sampledf_blankcip_T96$Treatment, bias.adjust = TRUE, type = "centroid")
anova(disp_blankcip_T96) # Sig p = 1.57e-05

stats_blankcip_T96 <- broom::tidy(p.adjust(permutest(betadisper(clr_euc_blankcip_T96, sampledf_blankcip_T96$Treatment, type = "centroid", 
                              bias.adjust = TRUE), pairwise=TRUE)$pairwise$permuted, method = 'fdr')) # Significant between blank-cip low and cip low-cip high

stats_blankcip_T96$Time <- 96

df_disp_blankcip_T96 <- data.frame(disp_blankcip_T96$distances)
colnames(df_disp_blankcip_T96)[colnames(df_disp_blankcip_T96) == 'disp_blankcip_T96.distances'] <- "Distances"
df_disp_blankcip_T96$Time <- 96

## Combine Blank Cip distance dataframes 
df_disp_blankcip <- rbind(df_disp_blankcip_T12, df_disp_blankcip_T24, df_disp_blankcip_T48, df_disp_blankcip_T96)
df_disp_blankcip$Group <- "Ciprofloxacin"

## Combine Blank Cip stats dataframes 
stats_blankcip <- rbind(stats_blankcip_T12, stats_blankcip_T24, stats_blankcip_T48, stats_blankcip_T96)

# Change column name from x to p.adj
colnames(stats_blankcip)[colnames(stats_blankcip) == 'x'] <- "p.adj"

# Add significance column
stats_blankcip <- add_significance(stats_blankcip, p.col = "p.adj")

# Add group column
stats_blankcip$Group <- "Ciprofloxacin"


## Combine all dispersion and stats dataframes ====
# Dispersion dataframe
df_disp <- rbind(df_disp_T0, df_disp_blankmix, df_disp_blankamp,df_disp_blankstrep, df_disp_blankcip)

rownames(df_disp) -> sample.names

samples <- as.character(sapply(strsplit(sample.names, "G7_T\\d+\\_"), `[`, 2))
samples <- as.character(sapply(strsplit(samples, "_[a-d]"), `[`, 1)) 

treatment <- as.character(sapply(strsplit(samples, "_"), `[`, 1)) 

## Add samples and treatment columns
df_disp$Samples <- samples
df_disp$Treatment <- treatment

## Change treatment column
df_disp$Treatment[df_disp$Treatment == "B2" |
                    df_disp$Treatment == "B3" |
                    df_disp$Treatment == "B4" |
                    df_disp$Treatment == "B5" |
                    df_disp$Treatment == "B6"] <- "Blank"

df_disp$Treatment[df_disp$Treatment == "M1"] <- "Mixture Low"
df_disp$Treatment[df_disp$Treatment == "M2"] <- "Mixture High"
df_disp$Treatment[df_disp$Treatment == "A1"] <- "Ampicillin Low"
df_disp$Treatment[df_disp$Treatment == "A2"] <- "Ampicillin High"
df_disp$Treatment[df_disp$Treatment == "S1"] <- "Streptomycin Low"
df_disp$Treatment[df_disp$Treatment == "S2"] <- "Streptomycin High"
df_disp$Treatment[df_disp$Treatment == "C1"] <- "Ciprofloxacin Low"
df_disp$Treatment[df_disp$Treatment == "C2"] <- "Ciprofloxacin High"

# Add Antibiotic column
df_disp$Antibiotic <- paste0(df_disp$Treatment)
df_disp$Antibiotic[df_disp$Antibiotic == "Mixture Low" |
                     df_disp$Antibiotic == "Mixture High"] <- "Mixture"
df_disp$Antibiotic[df_disp$Antibiotic == "Ampicillin Low" |
                     df_disp$Antibiotic == "Ampicillin High"] <- "Ampicillin"
df_disp$Antibiotic[df_disp$Antibiotic == "Streptomycin Low" |
                     df_disp$Antibiotic == "Streptomycin High"] <- "Streptomycin"
df_disp$Antibiotic[df_disp$Antibiotic == "Ciprofloxacin Low" |
                     df_disp$Antibiotic == "Ciprofloxacin High"] <- "Ciprofloxacin"


write.csv(df_disp_all, here::here("Data/08 - Beta Diversity - Output/06 - dispersion_dist.csv"))

# Statistics dataframe
df_disp_stats <- rbind(stats_T0, stats_blankmix, stats_blankamp, stats_blankstrep, stats_blankcip)

## Format for plotting 
### Add y.position
df_disp_stats$y.position <- 10

### Add group1 and group2 columns for stat_pvalue_manual (group 1 and 2 indicate x positions of bars)
df_disp_stats$group1 <- 1 # arbitrary for now
df_disp_stats$group2 <- 1

### Look into df_disp_stats and see which are significant
#### Blank-Mixture High T48
df_disp_stats[43,]$group1 <- 2.75
df_disp_stats[43,]$group2 <- 3.25

#### Ampicillin High-Ampicillin Low T48
df_disp_stats[55,]$group1 <- 3
df_disp_stats[55,]$group2 <- 3.25

#### Ampicillin High-Blank T48
df_disp_stats[56,]$group1 <- 2.75
df_disp_stats[56,]$group2 <- 3.25

#### Blank-Ciprofloxacin Low T24
df_disp_stats[77,]$group1 <- 1.75
df_disp_stats[77,]$group2 <- 2

#### Ciprofloxacin High-Ciprofloxacin Low T24
df_disp_stats[78,]$group1 <- 2
df_disp_stats[78,]$group2 <- 2.25

#### Blank-Ciprofloxacin Low T96
df_disp_stats[83,]$group1 <- 3.75
df_disp_stats[83,]$group2 <- 4

#### Ciprofloxacin High-Ciprofloxacin Low T96
df_disp_stats[84,]$group1 <- 4
df_disp_stats[84,]$group2 <- 4.25


## Plot ====
## Mixture
plot_disp_blankmix <-  ggplot(subset(df_disp, Group == "Mixture"), aes(x = as.factor(Time), y = Distances, color = Treatment))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("#666666", "#FF7856", "#FF3D3D"))
plot_disp_blankmix <- plot_disp_blankmix + geom_point(aes(color = Treatment), alpha = 0.5, 
                                    position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_blankmix <- plot_disp_blankmix + xlab("Time (Hours)") + ylab("Distance to Centroid") + ylim(1.75, 11)
plot_disp_blankmix <- plot_disp_blankmix + theme(axis.text = element_text(face = "bold", size = 11.5), 
                               axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  stat_pvalue_manual(subset(df_disp_stats, Group == "Mixture"), label = "p.adj.signif", hide.ns = TRUE)

## Ampicillin
plot_disp_blankamp <-  ggplot(subset(df_disp, Group == "Ampicillin"), aes(x = as.factor(Time), y = Distances, color = Treatment))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("#666666", "#9999FF" ,"#BC71BB"))
plot_disp_blankamp <- plot_disp_blankamp + geom_point(aes(color = Treatment), alpha = 0.5, 
                                                      position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_blankamp <- plot_disp_blankamp + xlab("Time (Hours)") + ylab("Distance to Centroid") + ylim(1.75, 11)
plot_disp_blankamp <- plot_disp_blankamp + theme(axis.text = element_text(face = "bold", size = 11.5), 
                                                 axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  stat_pvalue_manual(subset(df_disp_stats, Group == "Ampicillin"), label = "p.adj.signif", hide.ns = TRUE)

### Change y.position and run plot again
df_disp_stats[56,]$y.position <- 10.5

## Streptomycin
plot_disp_blankstrep <-  ggplot(subset(df_disp, Group == "Streptomycin"), aes(x = as.factor(Time), y = Distances, color = Treatment))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("#666666", "#3FA0FF", "#006DDB"))
plot_disp_blankstrep <- plot_disp_blankstrep + geom_point(aes(color = Treatment), alpha = 0.5, 
                                                      position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_blankstrep <- plot_disp_blankstrep + xlab("Time (Hours)") + ylab("Distance to Centroid") + ylim(1.75, 11)
plot_disp_blankstrep <- plot_disp_blankstrep + theme(axis.text = element_text(face = "bold", size = 11.5), 
                                                 axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  stat_pvalue_manual(subset(df_disp_stats, Group == "Streptomycin"), label = "p.adj.signif", hide.ns = TRUE)

## Ciprofloxacin
plot_disp_blankcip <-  ggplot(subset(df_disp, Group == "Ciprofloxacin"), aes(x = as.factor(Time), y = Distances, color = Treatment))  + 
  geom_boxplot(lwd = 1.25, outlier.colour = "NA") + 
  scale_color_manual(values = c("#666666", "#00BB00", "#008600"))
plot_disp_blankcip <- plot_disp_blankcip + geom_point(aes(color = Treatment), alpha = 0.5, 
                                                          position = position_jitterdodge(jitter.width = 0.1)) +theme_bw()
plot_disp_blankcip <- plot_disp_blankcip + xlab("Time (Hours)") + ylab("Distance to Centroid") + ylim(1.75, 11)
plot_disp_blankcip <- plot_disp_blankcip + theme(axis.text = element_text(face = "bold", size = 11.5), 
                                                     axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) + 
  guides(color = guide_legend(title = "Treatment")) + theme(strip.text = element_text(face = "bold", size = 12)) +
  stat_pvalue_manual(subset(df_disp_stats, Group == "Ciprofloxacin"), label = "p.adj.signif", hide.ns = TRUE)

### Change y.position and run plot again
df_disp_stats[78,]$y.position <- 10.5
df_disp_stats[84,]$y.position <- 10.5

## Save beta dispersion stats output
write.csv(df_disp_stats, here::here("Data/08 - Beta Diversity - Output/07 - dispersion_stats.csv"))


## Combine PCA and Dispersion plots ====
# Remove legend from PCA - treatment plot for easier plotting 
PCA_blankmix <- PCA_blankmix + theme(legend.position = "none") 
PCA_blankamp <- PCA_blankamp + theme(legend.position = "none") 
PCA_blankstrep <- PCA_blankstrep + theme(legend.position = "none") 
PCA_blankcip <- PCA_blankcip + theme(legend.position = "none")

grid.arrange(PCA_blankmix, PCA_blankamp, PCA_blankstrep, PCA_blankcip, ncol = 1) -> PCA_combined

grid.arrange(plot_disp_blankmix, plot_disp_blankamp, plot_disp_blankstrep, plot_disp_blankcip, ncol = 1) -> disp_combined

grid.arrange(PCA_combined, disp_combined, ncol = 2) -> plot_all

ggplot2::ggsave(here::here("Data/08 - Beta Diversity - Output/PCA_disp_final.png"), plot_all,
                height = 450, width = 650, units = "mm",
                scale = 0.5, dpi = 1000)
