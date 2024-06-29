# Author: Sunni Patton
# Date: 03/24/24 (last edited)
# Title: Differential Abundance Analysis 

## Set seed ====
set.seed(123)

## Load packages ====
library(dplyr)
library(ggplot2)
library(tidyr)
library(ANCOMBC)
library(phyloseq)
library(mia)
library(lme4)
library(lmerTest)

## Load data ====
readRDS(here::here("Data/09 - Differential Abundance - Output/ps.pruned.rds")) -> ps.pruned

## Compare T0 (pretreatment) to Antibiotic (no dose or time) ====

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.pruned)

## Set order for Treatment
tse$PretreatmentAbx <- factor(tse$PretreatmentAbx, levels = c("Pretreatment", "Blank", "Mixture", "Ampicillin",
                                                              "Streptomycin", "Ciprofloxacin"))

## ANCOMBC formula
### CoralID as random effect
output_abx <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "PretreatmentAbx", 
                           rand_formula = "(1 | CoralID)", p_adj_method = "fdr", group = "PretreatmentAbx", 
                           alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

saveRDS(output_abx,here::here("Data/09 - Differential Abundance - Output/01 - ANCOMBC2_output_abx.rds") )

## Save output 
res_prim <- output_abx$res

### Save specific output as dataframe 
#### Blank
df_blank <- data.frame(c(res_prim[1], res_prim[3], res_prim[27], res_prim[33]))
df_blank$Antibiotic <- "Blank"
colnames(df_blank)[2] <- "LFC"
colnames(df_blank)[3] <- "adj.p"
colnames(df_blank)[4] <- "Diff"

df_mix <- data.frame(c(res_prim[1], res_prim[4], res_prim[28], res_prim[34]))
df_mix$Antibiotic <- "Mixture"
colnames(df_mix)[2] <- "LFC"
colnames(df_mix)[3] <- "adj.p"
colnames(df_mix)[4] <- "Diff"

df_amp <- data.frame(c(res_prim[1], res_prim[5], res_prim[29], res_prim[35]))
df_amp$Antibiotic <- "Ampicillin"
colnames(df_amp)[2] <- "LFC"
colnames(df_amp)[3] <- "adj.p"
colnames(df_amp)[4] <- "Diff"

df_strep <- data.frame(c(res_prim[1], res_prim[6], res_prim[30], res_prim[36]))
df_strep$Antibiotic <- "Streptomycin"
colnames(df_strep)[2] <- "LFC"
colnames(df_strep)[3] <- "adj.p"
colnames(df_strep)[4] <- "Diff"

df_cip <- data.frame(c(res_prim[1], res_prim[7], res_prim[31], res_prim[37]))
df_cip$Antibiotic <- "Ciprofloxacin"
colnames(df_cip)[2] <- "LFC"
colnames(df_cip)[3] <- "adj.p"
colnames(df_cip)[4] <- "Diff"

### Merge dataframes
df_abx <- dplyr::bind_rows(df_blank, df_mix, df_amp, df_strep, df_cip)

### Add another taxon column for plot coloring
df_abx$Genus <- paste0(df_abx$taxon)
df_abx$Genus[df_abx$adj.p >= 0.05] <- "Other"
df_abx$Genus[df_abx$Genus == "ASV1"] <- "Campylobacterales Order (ASV1)" #E3211C
df_abx$Genus[df_abx$Genus == "ASV2"] <- "Alteromonas (ASV2)" #1F78B4
df_abx$Genus[df_abx$Genus == "ASV3"] <- "Helicobacteraceae Family (ASV3)" #33A02B
df_abx$Genus[df_abx$Genus == "ASV4"] <- "Phaeobacter (ASV4)"#FF7F00
df_abx$Genus[df_abx$Genus == "ASV5"] <- "P3OB-42 (ASV5)" #6A3D9A
df_abx$Genus[df_abx$Genus == "ASV9"] <- "[Caedibacter] taeniospiralis group (ASV9)" #FB9A99
df_abx$Genus[df_abx$Genus == "ASV12"] <- "Hyphomonadaceae Family (ASV12)" #B2DF8A
df_abx$Genus[df_abx$Genus == "ASV18"] <- "Alteromonadaceae Family (ASV18)" #FEC655

## Save dataframe
write.csv(df_abx, here::here("Data/09 - Differential Abundance - Output/01 - ANCOMBC2_output_abx.csv"))

## Set order for Taxa
df_abx$Genus <- factor(df_abx$Genus, levels = c("Campylobacterales Order (ASV1)", "Alteromonas (ASV2)",
                                                        "Helicobacteraceae Family (ASV3)", "Phaeobacter (ASV4)","P3OB-42 (ASV5)",
                                                        "[Caedibacter] taeniospiralis group (ASV9)","Hyphomonadaceae Family (ASV12)", 
                                                        "Alteromonadaceae Family (ASV18)", "Other"))
## Set order for Treatments
df_abx$Antibiotic <- factor(df_abx$Antibiotic, levels = c("Blank", "Mixture", "Ampicillin", "Streptomycin", "Ciprofloxacin"))

## Plot
vol_plot <- df_abx %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = Genus)) + 
  geom_point(size = 3.5, alpha = 0.8) + facet_wrap(~Antibiotic, ncol = 5) + scale_color_manual(values = c("#33A02B","#1F78B4","#E3211C",
                                                                                         "#FF7F00", "#FB9A99", "#CAB2D6", "#A9683F",
                                                                                         "#C89D6D", "darkgrey"))

vol_plot <- vol_plot + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
                     limits = c(-6, 6)) #+(0, 5)

vol_plot <- vol_plot + labs(title = "Antibiotic Treatments vs. Pretreatment (T0)") + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value")

ggplot2::ggsave(here::here("Data/09 - Differential Abundance - Output/01 - plot_PretreatAbx.png"), vol_plot,
                height = 250, width = 600, units = "mm",
                scale = 0.5, dpi = 1000)


## Compare T0 (pretreatment) to Treatment (no time) ====

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.pruned)

## Set order for Treatment
tse$Pretreatment <- factor(tse$Pretreatment, levels = c("Pretreatment", "Blank", "Mixture Low", "Mixture High",
                                                        "Ampicillin Low", "Ampicillin High", "Streptomycin Low", 
                                                        "Streptomycin High", "Ciprofloxacin Low", "Ciprofloxacin High"))

## ANCOMBC formula
### CoralID as random effect
output_treat <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "Pretreatment", 
                       rand_formula = "(1 | CoralID)", p_adj_method = "fdr", group = "Pretreatment",
                       alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

saveRDS(output_treat,here::here("Data/09 - Differential Abundance - Output/02 - ANCOMBC2_output_treat.rds") )

## Save output 
res_prim <- output_treat$res

### Save specific output as dataframe 
#### Blank
df_blank <- data.frame(c(res_prim[1], res_prim[3], res_prim[43], res_prim[53]))
df_blank$Treatment <- "Blank"
colnames(df_blank)[2] <- "LFC"
colnames(df_blank)[3] <- "adj.p"
colnames(df_blank)[4] <- "Diff"

#### Mixture Low
df_mixlow <- data.frame(c(res_prim[1], res_prim[4], res_prim[44], res_prim[54]))
df_mixlow$Treatment <- "Mixture Low"
colnames(df_mixlow)[2] <- "LFC"
colnames(df_mixlow)[3] <- "adj.p"
colnames(df_mixlow)[4] <- "Diff"

#### Mixture High
df_mixhigh <- data.frame(c(res_prim[1], res_prim[5], res_prim[45], res_prim[55]))
df_mixhigh$Treatment <- "Mixture High"
colnames(df_mixhigh)[2] <- "LFC"
colnames(df_mixhigh)[3] <- "adj.p"
colnames(df_mixhigh)[4] <- "Diff"

#### Ampicillin Low
df_amplow <- data.frame(c(res_prim[1], res_prim[6], res_prim[46], res_prim[56]))
df_amplow$Treatment <- "Ampicillin Low"
colnames(df_amplow)[2] <- "LFC"
colnames(df_amplow)[3] <- "adj.p"
colnames(df_amplow)[4] <- "Diff"

#### Ampicillin High
df_amphigh <- data.frame(c(res_prim[1], res_prim[7], res_prim[47], res_prim[57]))
df_amphigh$Treatment <- "Ampicillin High"
colnames(df_amphigh)[2] <- "LFC"
colnames(df_amphigh)[3] <- "adj.p"
colnames(df_amphigh)[4] <- "Diff"

#### Streptomycin Low
df_streplow <- data.frame(c(res_prim[1], res_prim[8], res_prim[48], res_prim[58]))
df_streplow$Treatment <- "Streptomycin Low"
colnames(df_streplow)[2] <- "LFC"
colnames(df_streplow)[3] <- "adj.p"
colnames(df_streplow)[4] <- "Diff"

#### Streptomycin High
df_strephigh <- data.frame(c(res_prim[1], res_prim[9], res_prim[49], res_prim[59]))
df_strephigh$Treatment <- "Streptomycin High"
colnames(df_strephigh)[2] <- "LFC"
colnames(df_strephigh)[3] <- "adj.p"
colnames(df_strephigh)[4] <- "Diff"

#### Ciprofloxacin Low
df_ciplow <- data.frame(c(res_prim[1], res_prim[10], res_prim[50], res_prim[60]))
df_ciplow$Treatment <- "Ciprofloxacin Low"
colnames(df_ciplow)[2] <- "LFC"
colnames(df_ciplow)[3] <- "adj.p"
colnames(df_ciplow)[4] <- "Diff"

#### Ciprofloxacin High
df_ciphigh <- data.frame(c(res_prim[1], res_prim[11], res_prim[51], res_prim[61]))
df_ciphigh$Treatment <- "Ciprofloxacin High"
colnames(df_ciphigh)[2] <- "LFC"
colnames(df_ciphigh)[3] <- "adj.p"
colnames(df_ciphigh)[4] <- "Diff"

### Merge dataframes
df_treat <- dplyr::bind_rows(df_blank, df_mixlow, df_mixhigh, df_amplow, df_amphigh, df_streplow, df_strephigh,
                             df_ciplow, df_ciphigh)

### Add another taxon column for plot coloring
df_treat$Genus <- paste0(df_treat$taxon)
df_treat$Genus[df_treat$adj.p >= 0.05] <- "Other"
df_treat$Genus[df_treat$Genus == "ASV1"] <- "Campylobacterales Order (ASV1)" #E3211C
df_treat$Genus[df_treat$Genus == "ASV2"] <- "Alteromonas (ASV2)" #1F78B4
df_treat$Genus[df_treat$Genus == "ASV3"] <- "Helicobacteraceae Family (ASV3)" #33A02B
df_treat$Genus[df_treat$Genus == "ASV4"] <- "Phaeobacter (ASV4)"#FF7F00
df_treat$Genus[df_treat$Genus == "ASV5"] <- "P3OB-42 (ASV5)" #6A3D9A
df_treat$Genus[df_treat$Genus == "ASV9"] <- "[Caedibacter] taeniospiralis group (ASV9)" #FB9A99
df_treat$Genus[df_treat$Genus == "ASV12"] <- "Hyphomonadaceae Family (ASV12)" #B2DF8A
df_treat$Genus[df_treat$Genus == "ASV14"] <- "Patescibacteria Phylum (ASV14)" #black 
df_treat$Genus[df_treat$Genus == "ASV18"] <- "Alteromonadaceae Family (ASV18)" #FEC655

## Save CSV
write.csv(df_treat, here::here("Data/09 - Differential Abundance - Output/02 - ANCOMBC2_output_treat.csv"))


## Set order for Taxa
df_treat$Genus <- factor(df_treat$Genus, levels = c("Campylobacterales Order (ASV1)", "Alteromonas (ASV2)",
                                                        "Helicobacteraceae Family (ASV3)", "Phaeobacter (ASV4)","P3OB-42 (ASV5)",
                                                        "[Caedibacter] taeniospiralis group (ASV9)","Hyphomonadaceae Family (ASV12)", 
                                                        "Patescibacteria Phylum (ASV14)", "Alteromonadaceae Family (ASV18)", "Other"))
## Set order for Treatments
df_treat$Treatment <- factor(df_treat$Treatment, levels = c("Blank", "Mixture Low", "Mixture High", "Ampicillin Low", "Streptomycin Low",
                                                            "Ciprofloxacin Low", "Ampicillin High", "Streptomycin High",
                                                            "Ciprofloxacin High"))

## Main Figure 4
vol_plot <- df_treat %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = Genus)) + 
  geom_point(size = 3.5, alpha = 0.8) + facet_wrap(~Treatment, ncol = 3, nrow = 3) + scale_color_manual(values = c("#33A02B", "#1F78B4","#E3211C",
                                                                                         "#FF7F00", "#FB9A99", "#CAB2D6", "#A9683F",
                                                                                         "black","#CFA77A", "darkgrey"))

vol_plot <- vol_plot + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
                     limits = c(-6, 6)) #+(0, 5)

vol_plot <- vol_plot + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value") + theme(strip.text = element_text(face = "bold", size = 12))

ggplot2::ggsave(here::here("Data/09 - Differential Abundance - Output/02.2 - plot_PretreatTreat.png"), vol_plot,
                height = 350, width = 550, units = "mm",
                scale = 0.5, dpi = 1000)


## Compare T0 (pretreatment) to PretreatmentTime ====

## Convert phyloseq object to tree summarized experiment
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps.pruned)

## Set order for Treatment
tse$PretreatmentTime <- factor(tse$PretreatmentTime, levels = c("Pretreatment", "Blank 12", "Blank 24", "Blank 48", "Blank 96",
                                                                "Mixture Low 12", "Mixture Low 24", "Mixture Low 48", "Mixture Low 96",
                                                                "Mixture High 12", "Mixture High 24", "Mixture High 48", "Mixture High 96", 
                                                                "Ampicillin Low 12", "Ampicillin Low 24", "Ampicillin Low 48", "Ampicillin Low 96",
                                                                "Ampicillin High 12", "Ampicillin High 24", "Ampicillin High 48", "Ampicillin High 96",
                                                                "Streptomycin Low 12","Streptomycin Low 24", "Streptomycin Low 48", "Streptomycin Low 96",
                                                                "Streptomycin High 12", "Streptomycin High 24", "Streptomycin High 48", "Streptomycin High 96",
                                                                "Ciprofloxacin Low 12", "Ciprofloxacin Low 24", "Ciprofloxacin Low 48", "Ciprofloxacin Low 96",
                                                                "Ciprofloxacin High 12", "Ciprofloxacin High 24", "Ciprofloxacin High 48", "Ciprofloxacin High 96"))

## ANCOMBC formula
### CoralID as random effect
output_treatTime <- ancombc2(data = tse, assay_name = "counts", tax_level = "ASV", fix_formula = "PretreatmentTime", 
                             rand_formula = "(1 | CoralID)", p_adj_method = "fdr", group = "PretreatmentTime",
                             alpha = 0.05, prv_cut = 0.1, neg_lb = FALSE) 

saveRDS(output_treatTime, here::here("Data/09 - Differential Abundance - Output/03 - ANCOMBC2_output_treatTime.rds"))

## Save output
res_prim <- output_treatTime$res

### Save specific output as dataframe 
#### Blank 
### T12
df_TreatTime_blankT12 <- data.frame(c(res_prim[1], res_prim[3], res_prim[151], res_prim[188]))
df_TreatTime_blankT12$Treatment <- "Blank"
df_TreatTime_blankT12$Time <- "12"
colnames(df_TreatTime_blankT12)[2] <- "LFC"
colnames(df_TreatTime_blankT12)[3] <- "adj.p"
colnames(df_TreatTime_blankT12)[4] <- "Diff"

### T24
df_TreatTime_blankT24 <- data.frame(c(res_prim[1], res_prim[4], res_prim[152], res_prim[189]))
df_TreatTime_blankT24$Treatment <- "Blank"
df_TreatTime_blankT24$Time <- "24"
colnames(df_TreatTime_blankT24)[2] <- "LFC"
colnames(df_TreatTime_blankT24)[3] <- "adj.p"
colnames(df_TreatTime_blankT24)[4] <- "Diff"


### T48
df_TreatTime_blankT48 <- data.frame(c(res_prim[1], res_prim[5], res_prim[153], res_prim[190]))
df_TreatTime_blankT48$Treatment <- "Blank"
df_TreatTime_blankT48$Time <- "48"
colnames(df_TreatTime_blankT48)[2] <- "LFC"
colnames(df_TreatTime_blankT48)[3] <- "adj.p"
colnames(df_TreatTime_blankT48)[4] <- "Diff"

### T96
df_TreatTime_blankT96 <- data.frame(c(res_prim[1], res_prim[6], res_prim[154], res_prim[191]))
df_TreatTime_blankT96$Treatment <- "Blank"
df_TreatTime_blankT96$Time <- "96"
colnames(df_TreatTime_blankT96)[2] <- "LFC"
colnames(df_TreatTime_blankT96)[3] <- "adj.p"
colnames(df_TreatTime_blankT96)[4] <- "Diff"


### Merge dataframes
df_TreatTime_blank <- dplyr::bind_rows(df_TreatTime_blankT12, df_TreatTime_blankT24, df_TreatTime_blankT48, df_TreatTime_blankT96)

#### Mixture Low
### T12
df_TreatTime_mixlowT12 <- data.frame(c(res_prim[1], res_prim[7], res_prim[155], res_prim[192]))
df_TreatTime_mixlowT12$Treatment <- "Mixture Low"
df_TreatTime_mixlowT12$Time <- "12"
colnames(df_TreatTime_mixlowT12)[2] <- "LFC"
colnames(df_TreatTime_mixlowT12)[3] <- "adj.p"
colnames(df_TreatTime_mixlowT12)[4] <- "Diff"

### T24
df_TreatTime_mixlowT24 <- data.frame(c(res_prim[1], res_prim[8], res_prim[156], res_prim[193]))
df_TreatTime_mixlowT24$Treatment <- "Mixture Low"
df_TreatTime_mixlowT24$Time <- "24"
colnames(df_TreatTime_mixlowT24)[2] <- "LFC"
colnames(df_TreatTime_mixlowT24)[3] <- "adj.p"
colnames(df_TreatTime_mixlowT24)[4] <- "Diff"

### T48
df_TreatTime_mixlowT48 <- data.frame(c(res_prim[1], res_prim[9], res_prim[157], res_prim[194]))
df_TreatTime_mixlowT48$Treatment <- "Mixture Low"
df_TreatTime_mixlowT48$Time <- "48"
colnames(df_TreatTime_mixlowT48)[2] <- "LFC"
colnames(df_TreatTime_mixlowT48)[3] <- "adj.p"
colnames(df_TreatTime_mixlowT48)[4] <- "Diff"

### T96
df_TreatTime_mixlowT96 <- data.frame(c(res_prim[1], res_prim[10], res_prim[158], res_prim[195]))
df_TreatTime_mixlowT96$Treatment <- "Mixture Low"
df_TreatTime_mixlowT96$Time <- "96"
colnames(df_TreatTime_mixlowT96)[2] <- "LFC"
colnames(df_TreatTime_mixlowT96)[3] <- "adj.p"
colnames(df_TreatTime_mixlowT96)[4] <- "Diff"

#### Mixture high
### T12
df_TreatTime_mixhighT12 <- data.frame(c(res_prim[1], res_prim[11], res_prim[159], res_prim[196]))
df_TreatTime_mixhighT12$Treatment <- "Mixture High"
df_TreatTime_mixhighT12$Time <- "12"
colnames(df_TreatTime_mixhighT12)[2] <- "LFC"
colnames(df_TreatTime_mixhighT12)[3] <- "adj.p"
colnames(df_TreatTime_mixhighT12)[4] <- "Diff"

### T24
df_TreatTime_mixhighT24 <- data.frame(c(res_prim[1], res_prim[12], res_prim[160], res_prim[197]))
df_TreatTime_mixhighT24$Treatment <- "Mixture High"
df_TreatTime_mixhighT24$Time <- "24"
colnames(df_TreatTime_mixhighT24)[2] <- "LFC"
colnames(df_TreatTime_mixhighT24)[3] <- "adj.p"
colnames(df_TreatTime_mixhighT24)[4] <- "Diff"

### T48
df_TreatTime_mixhighT48 <- data.frame(c(res_prim[1], res_prim[13], res_prim[161], res_prim[198]))
df_TreatTime_mixhighT48$Treatment <- "Mixture High"
df_TreatTime_mixhighT48$Time <- "48"
colnames(df_TreatTime_mixhighT48)[2] <- "LFC"
colnames(df_TreatTime_mixhighT48)[3] <- "adj.p"
colnames(df_TreatTime_mixhighT48)[4] <- "Diff"

### T96
df_TreatTime_mixhighT96 <- data.frame(c(res_prim[1], res_prim[14], res_prim[162], res_prim[199]))
df_TreatTime_mixhighT96$Treatment <- "Mixture High"
df_TreatTime_mixhighT96$Time <- "96"
colnames(df_TreatTime_mixhighT96)[2] <- "LFC"
colnames(df_TreatTime_mixhighT96)[3] <- "adj.p"
colnames(df_TreatTime_mixhighT96)[4] <- "Diff"

### Merge dataframes
df_TreatTime_mix <- dplyr::bind_rows(df_TreatTime_mixlowT12, df_TreatTime_mixlowT24, df_TreatTime_mixlowT48, df_TreatTime_mixlowT96,
                                     df_TreatTime_mixhighT12, df_TreatTime_mixhighT24, df_TreatTime_mixhighT48, df_TreatTime_mixhighT96)

#### Ampicillin low
### T12
df_TreatTime_amplowT12 <- data.frame(c(res_prim[1], res_prim[15], res_prim[163], res_prim[200]))
df_TreatTime_amplowT12$Treatment <- "Ampicillin Low"
df_TreatTime_amplowT12$Time <- "12"
colnames(df_TreatTime_amplowT12)[2] <- "LFC"
colnames(df_TreatTime_amplowT12)[3] <- "adj.p"
colnames(df_TreatTime_amplowT12)[4] <- "Diff"

### T24
df_TreatTime_amplowT24 <- data.frame(c(res_prim[1], res_prim[16], res_prim[164], res_prim[201]))
df_TreatTime_amplowT24$Treatment <- "Ampicillin Low"
df_TreatTime_amplowT24$Time <- "24"
colnames(df_TreatTime_amplowT24)[2] <- "LFC"
colnames(df_TreatTime_amplowT24)[3] <- "adj.p"
colnames(df_TreatTime_amplowT24)[4] <- "Diff"

### T48
df_TreatTime_amplowT48 <- data.frame(c(res_prim[1], res_prim[17], res_prim[165], res_prim[202]))
df_TreatTime_amplowT48$Treatment <- "Ampicillin Low"
df_TreatTime_amplowT48$Time <- "48"
colnames(df_TreatTime_amplowT48)[2] <- "LFC"
colnames(df_TreatTime_amplowT48)[3] <- "adj.p"
colnames(df_TreatTime_amplowT48)[4] <- "Diff"

### T96
df_TreatTime_amplowT96 <- data.frame(c(res_prim[1], res_prim[18], res_prim[166], res_prim[203]))
df_TreatTime_amplowT96$Treatment <- "Ampicillin Low"
df_TreatTime_amplowT96$Time <- "96"
colnames(df_TreatTime_amplowT96)[2] <- "LFC"
colnames(df_TreatTime_amplowT96)[3] <- "adj.p"
colnames(df_TreatTime_amplowT96)[4] <- "Diff"

#### Amp High T12
df_TreatTime_ampHighT12 <- data.frame(c(res_prim[1], res_prim[19], res_prim[167], res_prim[204]))
df_TreatTime_ampHighT12$Treatment <- "Ampicillin High"
df_TreatTime_ampHighT12$Time <- "12"
colnames(df_TreatTime_ampHighT12)[2] <- "LFC"
colnames(df_TreatTime_ampHighT12)[3] <- "adj.p"
colnames(df_TreatTime_ampHighT12)[4] <- "Diff"

#### Amp High T24
df_TreatTime_ampHighT24 <- data.frame(c(res_prim[1], res_prim[20], res_prim[168], res_prim[205]))
df_TreatTime_ampHighT24$Treatment <- "Ampicillin High"
df_TreatTime_ampHighT24$Time <- "24"
colnames(df_TreatTime_ampHighT24)[2] <- "LFC"
colnames(df_TreatTime_ampHighT24)[3] <- "adj.p"
colnames(df_TreatTime_ampHighT24)[4] <- "Diff"

#### Amp High T48
df_TreatTime_ampHighT48 <- data.frame(c(res_prim[1], res_prim[21], res_prim[169], res_prim[206]))
df_TreatTime_ampHighT48$Treatment <- "Ampicillin High"
df_TreatTime_ampHighT48$Time <- "48"
colnames(df_TreatTime_ampHighT48)[2] <- "LFC"
colnames(df_TreatTime_ampHighT48)[3] <- "adj.p"
colnames(df_TreatTime_ampHighT48)[4] <- "Diff"

#### Amp High T96
df_TreatTime_ampHighT96 <- data.frame(c(res_prim[1], res_prim[22], res_prim[170], res_prim[207]))
df_TreatTime_ampHighT96$Treatment <- "Ampicillin High"
df_TreatTime_ampHighT96$Time <- "96"
colnames(df_TreatTime_ampHighT96)[2] <- "LFC"
colnames(df_TreatTime_ampHighT96)[3] <- "adj.p"
colnames(df_TreatTime_ampHighT96)[4] <- "Diff"

### Merge dataframes
df_TreatTime_amp <- dplyr::bind_rows(df_TreatTime_amplowT12, df_TreatTime_amplowT24, df_TreatTime_amplowT48, df_TreatTime_amplowT96,
                                     df_TreatTime_ampHighT12, df_TreatTime_ampHighT24, df_TreatTime_ampHighT48, df_TreatTime_ampHighT96)

#### Strep Low T12
df_TreatTime_strepLowT12 <- data.frame(c(res_prim[1], res_prim[23], res_prim[171], res_prim[208]))
df_TreatTime_strepLowT12$Treatment <- "Streptomycin Low"
df_TreatTime_strepLowT12$Time <- "12"
colnames(df_TreatTime_strepLowT12)[2] <- "LFC"
colnames(df_TreatTime_strepLowT12)[3] <- "adj.p"
colnames(df_TreatTime_strepLowT12)[4] <- "Diff"

#### Strep Low T24
df_TreatTime_strepLowT24 <- data.frame(c(res_prim[1], res_prim[24], res_prim[172], res_prim[209]))
df_TreatTime_strepLowT24$Treatment <- "Streptomycin Low"
df_TreatTime_strepLowT24$Time <- "24"
colnames(df_TreatTime_strepLowT24)[2] <- "LFC"
colnames(df_TreatTime_strepLowT24)[3] <- "adj.p"
colnames(df_TreatTime_strepLowT24)[4] <- "Diff"

#### Strep Low T48
df_TreatTime_strepLowT48 <- data.frame(c(res_prim[1], res_prim[25], res_prim[173], res_prim[210]))
df_TreatTime_strepLowT48$Treatment <- "Streptomycin Low"
df_TreatTime_strepLowT48$Time <- "48"
colnames(df_TreatTime_strepLowT48)[2] <- "LFC"
colnames(df_TreatTime_strepLowT48)[3] <- "adj.p"
colnames(df_TreatTime_strepLowT48)[4] <- "Diff"

#### Strep Low T96
df_TreatTime_strepLowT96 <- data.frame(c(res_prim[1], res_prim[26], res_prim[174], res_prim[211]))
df_TreatTime_strepLowT96$Treatment <- "Streptomycin Low"
df_TreatTime_strepLowT96$Time <- "96"
colnames(df_TreatTime_strepLowT96)[2] <- "LFC"
colnames(df_TreatTime_strepLowT96)[3] <- "adj.p"
colnames(df_TreatTime_strepLowT96)[4] <- "Diff"

#### Strep High T12
df_TreatTime_strepHighT12 <- data.frame(c(res_prim[1], res_prim[27], res_prim[175], res_prim[212]))
df_TreatTime_strepHighT12$Treatment <- "Streptomycin High"
df_TreatTime_strepHighT12$Time <- "12"
colnames(df_TreatTime_strepHighT12)[2] <- "LFC"
colnames(df_TreatTime_strepHighT12)[3] <- "adj.p"
colnames(df_TreatTime_strepHighT12)[4] <- "Diff"

#### Strep High T24
df_TreatTime_strepHighT24 <- data.frame(c(res_prim[1], res_prim[28], res_prim[176], res_prim[213]))
df_TreatTime_strepHighT24$Treatment <- "Streptomycin High"
df_TreatTime_strepHighT24$Time <- "24"
colnames(df_TreatTime_strepHighT24)[2] <- "LFC"
colnames(df_TreatTime_strepHighT24)[3] <- "adj.p"
colnames(df_TreatTime_strepHighT24)[4] <- "Diff"

#### Strep High T48
df_TreatTime_strepHighT48 <- data.frame(c(res_prim[1], res_prim[29], res_prim[177], res_prim[214]))
df_TreatTime_strepHighT48$Treatment <- "Streptomycin High"
df_TreatTime_strepHighT48$Time <- "48"
colnames(df_TreatTime_strepHighT48)[2] <- "LFC"
colnames(df_TreatTime_strepHighT48)[3] <- "adj.p"
colnames(df_TreatTime_strepHighT48)[4] <- "Diff"

#### Strep High T96
df_TreatTime_strepHighT96 <- data.frame(c(res_prim[1], res_prim[30], res_prim[178], res_prim[215]))
df_TreatTime_strepHighT96$Treatment <- "Streptomycin High"
df_TreatTime_strepHighT96$Time <- "96"
colnames(df_TreatTime_strepHighT96)[2] <- "LFC"
colnames(df_TreatTime_strepHighT96)[3] <- "adj.p"
colnames(df_TreatTime_strepHighT96)[4] <- "Diff"

### Merge dataframes
df_TreatTime_strep <- dplyr::bind_rows(df_TreatTime_strepLowT12, df_TreatTime_strepLowT24, df_TreatTime_strepLowT48, df_TreatTime_strepLowT96,
                                       df_TreatTime_strepHighT12, df_TreatTime_strepHighT24, df_TreatTime_strepHighT48, df_TreatTime_strepHighT96)

#### Cip Low T12
df_TreatTime_cipLowT12 <- data.frame(c(res_prim[1], res_prim[31], res_prim[179], res_prim[216]))
df_TreatTime_cipLowT12$Treatment <- "Ciprofloxacin Low"
df_TreatTime_cipLowT12$Time <- "12"
colnames(df_TreatTime_cipLowT12)[2] <- "LFC"
colnames(df_TreatTime_cipLowT12)[3] <- "adj.p"
colnames(df_TreatTime_cipLowT12)[4] <- "Diff"

#### Cip Low T24
df_TreatTime_cipLowT24 <- data.frame(c(res_prim[1], res_prim[32], res_prim[180], res_prim[217]))
df_TreatTime_cipLowT24$Treatment <- "Ciprofloxacin Low"
df_TreatTime_cipLowT24$Time <- "24"
colnames(df_TreatTime_cipLowT24)[2] <- "LFC"
colnames(df_TreatTime_cipLowT24)[3] <- "adj.p"
colnames(df_TreatTime_cipLowT24)[4] <- "Diff"

#### Cip Low T48
df_TreatTime_cipLowT48 <- data.frame(c(res_prim[1], res_prim[33], res_prim[181], res_prim[218]))
df_TreatTime_cipLowT48$Treatment <- "Ciprofloxacin Low"
df_TreatTime_cipLowT48$Time <- "48"
colnames(df_TreatTime_cipLowT48)[2] <- "LFC"
colnames(df_TreatTime_cipLowT48)[3] <- "adj.p"
colnames(df_TreatTime_cipLowT48)[4] <- "Diff"

#### Cip Low T96
df_TreatTime_cipLowT96 <- data.frame(c(res_prim[1], res_prim[34], res_prim[182], res_prim[219]))
df_TreatTime_cipLowT96$Treatment <- "Ciprofloxacin Low"
df_TreatTime_cipLowT96$Time <- "96"
colnames(df_TreatTime_cipLowT96)[2] <- "LFC"
colnames(df_TreatTime_cipLowT96)[3] <- "adj.p"
colnames(df_TreatTime_cipLowT96)[4] <- "Diff"

#### Cip High T12
df_TreatTime_cipHighT12 <- data.frame(c(res_prim[1], res_prim[35], res_prim[183], res_prim[220]))
df_TreatTime_cipHighT12$Treatment <- "Ciprofloxacin High"
df_TreatTime_cipHighT12$Time <- "12"
colnames(df_TreatTime_cipHighT12)[2] <- "LFC"
colnames(df_TreatTime_cipHighT12)[3] <- "adj.p"
colnames(df_TreatTime_cipHighT12)[4] <- "Diff"

#### Cip High T24
df_TreatTime_cipHighT24 <- data.frame(c(res_prim[1], res_prim[36], res_prim[184], res_prim[221]))
df_TreatTime_cipHighT24$Treatment <- "Ciprofloxacin High"
df_TreatTime_cipHighT24$Time <- "24"
colnames(df_TreatTime_cipHighT24)[2] <- "LFC"
colnames(df_TreatTime_cipHighT24)[3] <- "adj.p"
colnames(df_TreatTime_cipHighT24)[4] <- "Diff"

#### Cip High T48
df_TreatTime_cipHighT48 <- data.frame(c(res_prim[1], res_prim[37], res_prim[185], res_prim[222]))
df_TreatTime_cipHighT48$Treatment <- "Ciprofloxacin High"
df_TreatTime_cipHighT48$Time <- "48"
colnames(df_TreatTime_cipHighT48)[2] <- "LFC"
colnames(df_TreatTime_cipHighT48)[3] <- "adj.p"
colnames(df_TreatTime_cipHighT48)[4] <- "Diff"

#### Cip High T96
df_TreatTime_cipHighT96 <- data.frame(c(res_prim[1], res_prim[38], res_prim[186], res_prim[223]))
df_TreatTime_cipHighT96$Treatment <- "Ciprofloxacin High"
df_TreatTime_cipHighT96$Time <- "96"
colnames(df_TreatTime_cipHighT96)[2] <- "LFC"
colnames(df_TreatTime_cipHighT96)[3] <- "adj.p"
colnames(df_TreatTime_cipHighT96)[4] <- "Diff"

### Merge dataframes
df_TreatTime_cip <- dplyr::bind_rows(df_TreatTime_cipLowT12, df_TreatTime_cipLowT24, df_TreatTime_cipLowT48, df_TreatTime_cipLowT96,
                                     df_TreatTime_cipHighT12, df_TreatTime_cipHighT24, df_TreatTime_cipHighT48, df_TreatTime_cipHighT96)


### Merge dataframes
df_TreatTime <- dplyr::bind_rows(df_TreatTime_blank, df_TreatTime_mix, df_TreatTime_amp, df_TreatTime_strep, df_TreatTime_cip)

### Add another taxon column for plot coloring
df_TreatTime$Genus <- paste0(df_TreatTime$taxon)
df_TreatTime$Genus[df_TreatTime$adj.p >= 0.05] <- "Other"
df_TreatTime$Genus[df_TreatTime$Genus == "ASV1"] <- "Campylobacterales Order (ASV1)" #E3211C
df_TreatTime$Genus[df_TreatTime$Genus == "ASV2"] <- "Alteromonas (ASV2)" #1F78B4
df_TreatTime$Genus[df_TreatTime$Genus == "ASV3"] <- "Helicobacteraceae Family (ASV3)" #33A02B
df_TreatTime$Genus[df_TreatTime$Genus == "ASV4"] <- "Nautella (ASV4)"#FF7F00
df_TreatTime$Genus[df_TreatTime$Genus == "ASV5"] <- "P3OB-42 (ASV5)" #6A3D9A
df_TreatTime$Genus[df_TreatTime$Genus == "ASV9"] <- "[Caedibacter] taeniospiralis group (ASV9)" #FB9A99
df_TreatTime$Genus[df_TreatTime$Genus == "ASV10"] <- "[Caedibacter] taeniospiralis group (ASV10)" #F1B6DA
df_TreatTime$Genus[df_TreatTime$Genus == "ASV11"] <- "Pseudoalteromonas (ASV11)" #0E6749
df_TreatTime$Genus[df_TreatTime$Genus == "ASV12"] <- "Hyphomonadaceae Family (ASV12)" #B2DF8A
df_TreatTime$Genus[df_TreatTime$Genus == "ASV14"] <- "Patescibacteria (ASV14)" #black
df_TreatTime$Genus[df_TreatTime$Genus == "ASV18"] <- "Alteromonadaceae Family (ASV18)" #FEC655
df_TreatTime$Genus[df_TreatTime$Genus == "ASV23"] <- "Tetragenococcus (ASV23)" #D08743
df_TreatTime$Genus[df_TreatTime$Genus == "ASV38"] <- "Idiomarina (ASV38)" #FF4757
df_TreatTime$Genus[df_TreatTime$Genus == "ASV39"] <- "Staphylococcus (ASV39)" #97D9D5



### Save dataframe
write.csv(df_TreatTime, here::here("Data/09 - Differential Abundance - Output/03 - ANCOMBC2_output_treatTime.csv"))


## Set order for Taxa
df_TreatTime$Genus <- factor(df_TreatTime$Genus, levels = c("Campylobacterales Order (ASV1)", "Alteromonas (ASV2)",
                                                            "Helicobacteraceae Family (ASV3)", "Nautella (ASV4)","P3OB-42 (ASV5)",
                                                            "[Caedibacter] taeniospiralis group (ASV9)", "[Caedibacter] taeniospiralis group (ASV10)",
                                                            "Pseudoalteromonas (ASV11)", "Hyphomonadaceae Family (ASV12)", "Patescibacteria (ASV14)", 
                                                            "Alteromonadaceae Family (ASV18)", "Tetragenococcus (ASV23)", "Idiomarina (ASV38)",
                                                            "Staphylococcus (ASV39)", "Other"))
## Set order for Treatments
df_TreatTime$Treatment <- factor(df_TreatTime$Treatment, levels = c("Blank", "Mixture Low", "Mixture High", "Ampicillin Low", "Ampicillin High", "Streptomycin Low",
                                                                    "Streptomycin High", "Ciprofloxacin Low", "Ciprofloxacin High"))

vol_plot <- df_TreatTime %>%
  ggplot(aes(x = LFC,
             y = -log10(adj.p),
             color = Genus)) + 
  geom_point(size = 3.5, alpha = 0.8) + facet_grid(cols = vars(Time), rows = vars(Treatment)) + scale_color_manual(values = c("#E3211C", "#1F78B4","#33A02B",
                                                                                                                              "#FF7F00", "#6A3D9A", "#FB9A99", "#F1B6DA", 
                                                                                                                              "#0E6749", "#B2DF8A", "black", "#FEC655", "#D08743",
                                                                                                                              "#FF4757", "#97D9D5", "darkgrey", "cyan", "hotpink"))
## Supplemental Figure 5
vol_plot <- vol_plot + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed") 
#scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
#                  limits = c(-6, 6)) #+(0, 5)

vol_plot <- vol_plot + labs(title = "Antibiotic Treatments vs. Pretreatment (T0)") + theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 11.5), 
        axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) +
  xlab("Log Fold Change") + ylab("-log10 Adjusted p value")

ggplot2::ggsave(here::here("Data/09 - Differential Abundance - Output/03 - plot_PretreatTime.png"), vol_plot,
                height = 600, width = 500, units = "mm",
                scale = 0.5, dpi = 1000)
