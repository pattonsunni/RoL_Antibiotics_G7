# Author: Sunni Patton
# Last edited: 6/3/24 
# Title: Relative Abundance

## Set seed ====
set.seed(123)

## Load libraries ====
library(microViz)
library(dplyr)
library(phyloseq)
library(ggpubr)

## Load data ====
readRDS(here::here("Data/05 - Phyloseq - Output/ps.rare.rds")) -> ps.rare

## Validate phyloseq object ====
tax_fix(
  ps.rare,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps.rare

ps.rare <- phyloseq_validate(ps.rare) # NAs detected

## Change taxa names ====
# Agree with current literature
taxa_change <- data.frame(tax_table(ps.rare))

for (i in 1:nrow(taxa_change)){
  
  if (taxa_change[i,6] == "MD3-55"){
    taxa_change[i, 6:7] <- "Aquarickettsia"
  }
  if (taxa_change[i,6] == "Nautella"){
    taxa_change[i, 6] <- "Phaeobacter"
  }
  if (taxa_change[i, 7] == "Nautella Genus"){
    taxa_change[i, 7] <- "Phaeobacter Genus"
  } 
}

tax_table(ps.rare) <- as.matrix(taxa_change)


## Remove B1 samples from analysis ====
ps.rare <- subset_samples(ps.rare, TreatmentTank != "B1")


## Relative abundance ps object ====
ps.rare.trans <- transform_sample_counts(ps.rare, function(OTU) OTU/sum(OTU))

## Identify top 10 most abundant taxa
top10 <- names(sort(taxa_sums(ps.rare.trans), decreasing = TRUE))[1:10] 

## Prune top 10
ps.rare.top10 <- prune_taxa(top10, ps.rare.trans)

## Relative abundance ====

## Supplemental Figure 4
plot.relAbund.all <- plot_bar(ps.rare.top10, x="CoralID", fill="Genus") + 
  facet_grid(vars(Time), vars(Antibiotic), scales = "free_x") 

plot.relAbund.all <- plot.relAbund.all + theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 11.5), 
        axis.text.y = element_text(face = "bold", size = 11.5), title = element_text(face = "bold")) + xlab("Samples") + labs(title = "Relative Abundance")


ggplot2::ggsave(here::here("Data/06 - Relative Abundance - Output/01 - plot_relabund_all.png"), plot.relAbund.all,
                height = 450, width = 900, units = "mm",
                scale = 0.5, dpi = 1000)

## Create dataframe for mean relative abundance Mixture Low ====
### Subset by treatment 
ps.mix.low <- subset_samples(ps.rare.top10, (Treatment == "Mixture Low"))

### Subset by time 
ps.mix.low.T0 <- subset_samples(ps.mix.low, (Time == "0"))
ps.mix.low.T12 <- subset_samples(ps.mix.low, (Time == "12"))
ps.mix.low.T24 <- subset_samples(ps.mix.low, (Time == "24"))
ps.mix.low.T48 <- subset_samples(ps.mix.low, (Time == "48"))
ps.mix.low.T96 <- subset_samples(ps.mix.low, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
mix.low.T0.df <- as.data.frame(ps.mix.low.T0@otu_table)
mix.low.T12.df <- as.data.frame(ps.mix.low.T12@otu_table)
mix.low.T24.df <- as.data.frame(ps.mix.low.T24@otu_table)
mix.low.T48.df <- as.data.frame(ps.mix.low.T48@otu_table)
mix.low.T96.df <- as.data.frame(ps.mix.low.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(mix.low.T0.df)

mix.low.T0.MRA <- data.frame(ASV = taxaOrder,
                             MeanRelAbund = c(mean(mix.low.T0.df$ASV4)*100,
                                              mean(mix.low.T0.df$ASV8)*100,
                                              mean(mix.low.T0.df$ASV9)*100,
                                              mean(mix.low.T0.df$ASV10)*100,
                                              mean(mix.low.T0.df$ASV2)*100,
                                              mean(mix.low.T0.df$ASV6)*100,
                                              mean(mix.low.T0.df$ASV5)*100,
                                              mean(mix.low.T0.df$ASV1)*100,
                                              mean(mix.low.T0.df$ASV3)*100,
                                              mean(mix.low.T0.df$ASV7)*100), 
                             SD = c(sd(mix.low.T0.df$ASV4)*100,
                                    sd(mix.low.T0.df$ASV8)*100,
                                    sd(mix.low.T0.df$ASV9)*100,
                                    sd(mix.low.T0.df$ASV10)*100,
                                    sd(mix.low.T0.df$ASV2)*100,
                                    sd(mix.low.T0.df$ASV6)*100,
                                    sd(mix.low.T0.df$ASV5)*100,
                                    sd(mix.low.T0.df$ASV1)*100,
                                    sd(mix.low.T0.df$ASV3)*100,
                                    sd(mix.low.T0.df$ASV7)*100),
                             Treatment = "Mixture Low",
                             Time = "0",
                             Taxon = colnames(mix.low.T0.df),
                             SamplesCombined = nrow(mix.low.T0.df))

### T12 dataframe
taxaOrder <- colnames(mix.low.T12.df)

mix.low.T12.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(mix.low.T12.df$ASV4)*100,
                                               mean(mix.low.T12.df$ASV8)*100,
                                               mean(mix.low.T12.df$ASV9)*100,
                                               mean(mix.low.T12.df$ASV10)*100,
                                               mean(mix.low.T12.df$ASV2)*100,
                                               mean(mix.low.T12.df$ASV6)*100,
                                               mean(mix.low.T12.df$ASV5)*100,
                                               mean(mix.low.T12.df$ASV1)*100,
                                               mean(mix.low.T12.df$ASV3)*100,
                                               mean(mix.low.T12.df$ASV7)*100), 
                              SD = c(sd(mix.low.T12.df$ASV4)*100,
                                     sd(mix.low.T12.df$ASV8)*100,
                                     sd(mix.low.T12.df$ASV9)*100,
                                     sd(mix.low.T12.df$ASV10)*100,
                                     sd(mix.low.T12.df$ASV2)*100,
                                     sd(mix.low.T12.df$ASV6)*100,
                                     sd(mix.low.T12.df$ASV5)*100,
                                     sd(mix.low.T12.df$ASV1)*100,
                                     sd(mix.low.T12.df$ASV3)*100,
                                     sd(mix.low.T12.df$ASV7)*100),
                              Treatment = "Mixture Low",
                              Time = "12",
                              Taxon = colnames(mix.low.T12.df),
                              SamplesCombined = nrow(mix.low.T12.df))

# T24 dataframe
taxaOrder <- colnames(mix.low.T24.df)

mix.low.T24.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(mix.low.T24.df$ASV4)*100,
                                               mean(mix.low.T24.df$ASV8)*100,
                                               mean(mix.low.T24.df$ASV9)*100,
                                               mean(mix.low.T24.df$ASV10)*100,
                                               mean(mix.low.T24.df$ASV2)*100,
                                               mean(mix.low.T24.df$ASV6)*100,
                                               mean(mix.low.T24.df$ASV5)*100,
                                               mean(mix.low.T24.df$ASV1)*100,
                                               mean(mix.low.T24.df$ASV3)*100,
                                               mean(mix.low.T24.df$ASV7)*100), 
                              SD = c(sd(mix.low.T24.df$ASV4)*100,
                                     sd(mix.low.T24.df$ASV8)*100,
                                     sd(mix.low.T24.df$ASV9)*100,
                                     sd(mix.low.T24.df$ASV10)*100,
                                     sd(mix.low.T24.df$ASV2)*100,
                                     sd(mix.low.T24.df$ASV6)*100,
                                     sd(mix.low.T24.df$ASV5)*100,
                                     sd(mix.low.T24.df$ASV1)*100,
                                     sd(mix.low.T24.df$ASV3)*100,
                                     sd(mix.low.T24.df$ASV7)*100),
                              Treatment = "Mixture Low",
                              Time = "24",
                              Taxon = colnames(mix.low.T24.df),
                              SamplesCombined = nrow(mix.low.T24.df))

# T48 dataframe
taxaOrder <- colnames(mix.low.T48.df)

mix.low.T48.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(mix.low.T48.df$ASV4)*100,
                                               mean(mix.low.T48.df$ASV8)*100,
                                               mean(mix.low.T48.df$ASV9)*100,
                                               mean(mix.low.T48.df$ASV10)*100,
                                               mean(mix.low.T48.df$ASV2)*100,
                                               mean(mix.low.T48.df$ASV6)*100,
                                               mean(mix.low.T48.df$ASV5)*100,
                                               mean(mix.low.T48.df$ASV1)*100,
                                               mean(mix.low.T48.df$ASV3)*100,
                                               mean(mix.low.T48.df$ASV7)*100), 
                              SD = c(sd(mix.low.T48.df$ASV4)*100,
                                     sd(mix.low.T48.df$ASV8)*100,
                                     sd(mix.low.T48.df$ASV9)*100,
                                     sd(mix.low.T48.df$ASV10)*100,
                                     sd(mix.low.T48.df$ASV2)*100,
                                     sd(mix.low.T48.df$ASV6)*100,
                                     sd(mix.low.T48.df$ASV5)*100,
                                     sd(mix.low.T48.df$ASV1)*100,
                                     sd(mix.low.T48.df$ASV3)*100,
                                     sd(mix.low.T48.df$ASV7)*100),
                              Treatment = "Mixture Low",
                              Time = "48",
                              Taxon = colnames(mix.low.T48.df),
                              SamplesCombined = nrow(mix.low.T48.df))

# T96 dataframe
taxaOrder <- colnames(mix.low.T96.df)

mix.low.T96.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(mix.low.T96.df$ASV4)*100,
                                               mean(mix.low.T96.df$ASV8)*100,
                                               mean(mix.low.T96.df$ASV9)*100,
                                               mean(mix.low.T96.df$ASV10)*100,
                                               mean(mix.low.T96.df$ASV2)*100,
                                               mean(mix.low.T96.df$ASV6)*100,
                                               mean(mix.low.T96.df$ASV5)*100,
                                               mean(mix.low.T96.df$ASV1)*100,
                                               mean(mix.low.T96.df$ASV3)*100,
                                               mean(mix.low.T96.df$ASV7)*100), 
                              SD = c(sd(mix.low.T96.df$ASV4)*100,
                                     sd(mix.low.T96.df$ASV8)*100,
                                     sd(mix.low.T96.df$ASV9)*100,
                                     sd(mix.low.T96.df$ASV10)*100,
                                     sd(mix.low.T96.df$ASV2)*100,
                                     sd(mix.low.T96.df$ASV6)*100,
                                     sd(mix.low.T96.df$ASV5)*100,
                                     sd(mix.low.T96.df$ASV1)*100,
                                     sd(mix.low.T96.df$ASV3)*100,
                                     sd(mix.low.T96.df$ASV7)*100),
                              Treatment = "Mixture Low",
                              Time = "96",
                              Taxon = colnames(mix.low.T96.df),
                              SamplesCombined = nrow(mix.low.T96.df))

## Create dataframe for mean relative abundance Mixture High ====
### Subset by treatment 
ps.mix.high <- subset_samples(ps.rare.top10, (Treatment == "Mixture High"))

### Subset by time 
ps.mix.high.T0 <- subset_samples(ps.mix.high, (Time == "0"))
ps.mix.high.T12 <- subset_samples(ps.mix.high, (Time == "12"))
ps.mix.high.T24 <- subset_samples(ps.mix.high, (Time == "24"))
ps.mix.high.T48 <- subset_samples(ps.mix.high, (Time == "48"))
ps.mix.high.T96 <- subset_samples(ps.mix.high, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
mix.high.T0.df <- as.data.frame(ps.mix.high.T0@otu_table)
mix.high.T12.df <- as.data.frame(ps.mix.high.T12@otu_table)
mix.high.T24.df <- as.data.frame(ps.mix.high.T24@otu_table)
mix.high.T48.df <- as.data.frame(ps.mix.high.T48@otu_table)
mix.high.T96.df <- as.data.frame(ps.mix.high.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(mix.high.T0.df)

mix.high.T0.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(mix.high.T0.df$ASV4)*100,
                                               mean(mix.high.T0.df$ASV8)*100,
                                               mean(mix.high.T0.df$ASV9)*100,
                                               mean(mix.high.T0.df$ASV10)*100,
                                               mean(mix.high.T0.df$ASV2)*100,
                                               mean(mix.high.T0.df$ASV6)*100,
                                               mean(mix.high.T0.df$ASV5)*100,
                                               mean(mix.high.T0.df$ASV1)*100,
                                               mean(mix.high.T0.df$ASV3)*100,
                                               mean(mix.high.T0.df$ASV7)*100), 
                              SD = c(sd(mix.high.T0.df$ASV4)*100,
                                     sd(mix.high.T0.df$ASV8)*100,
                                     sd(mix.high.T0.df$ASV9)*100,
                                     sd(mix.high.T0.df$ASV10)*100,
                                     sd(mix.high.T0.df$ASV2)*100,
                                     sd(mix.high.T0.df$ASV6)*100,
                                     sd(mix.high.T0.df$ASV5)*100,
                                     sd(mix.high.T0.df$ASV1)*100,
                                     sd(mix.high.T0.df$ASV3)*100,
                                     sd(mix.high.T0.df$ASV7)*100),
                              Treatment = "Mixture High",
                              Time = "0",
                              Taxon = colnames(mix.high.T0.df),
                              SamplesCombined = nrow(mix.high.T0.df))

### T12 dataframe
taxaOrder <- colnames(mix.high.T12.df)

mix.high.T12.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(mix.high.T12.df$ASV4)*100,
                                                mean(mix.high.T12.df$ASV8)*100,
                                                mean(mix.high.T12.df$ASV9)*100,
                                                mean(mix.high.T12.df$ASV10)*100,
                                                mean(mix.high.T12.df$ASV2)*100,
                                                mean(mix.high.T12.df$ASV6)*100,
                                                mean(mix.high.T12.df$ASV5)*100,
                                                mean(mix.high.T12.df$ASV1)*100,
                                                mean(mix.high.T12.df$ASV3)*100,
                                                mean(mix.high.T12.df$ASV7)*100), 
                               SD = c(sd(mix.high.T12.df$ASV4)*100,
                                      sd(mix.high.T12.df$ASV8)*100,
                                      sd(mix.high.T12.df$ASV9)*100,
                                      sd(mix.high.T12.df$ASV10)*100,
                                      sd(mix.high.T12.df$ASV2)*100,
                                      sd(mix.high.T12.df$ASV6)*100,
                                      sd(mix.high.T12.df$ASV5)*100,
                                      sd(mix.high.T12.df$ASV1)*100,
                                      sd(mix.high.T12.df$ASV3)*100,
                                      sd(mix.high.T12.df$ASV7)*100),
                               Treatment = "Mixture High",
                               Time = "12",
                               Taxon = colnames(mix.high.T12.df),
                               SamplesCombined = nrow(mix.high.T12.df))

# T24 dataframe
taxaOrder <- colnames(mix.high.T24.df)

mix.high.T24.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(mix.high.T24.df$ASV4)*100,
                                                mean(mix.high.T24.df$ASV8)*100,
                                                mean(mix.high.T24.df$ASV9)*100,
                                                mean(mix.high.T24.df$ASV10)*100,
                                                mean(mix.high.T24.df$ASV2)*100,
                                                mean(mix.high.T24.df$ASV6)*100,
                                                mean(mix.high.T24.df$ASV5)*100,
                                                mean(mix.high.T24.df$ASV1)*100,
                                                mean(mix.high.T24.df$ASV3)*100,
                                                mean(mix.high.T24.df$ASV7)*100), 
                               SD = c(sd(mix.high.T24.df$ASV4)*100,
                                      sd(mix.high.T24.df$ASV8)*100,
                                      sd(mix.high.T24.df$ASV9)*100,
                                      sd(mix.high.T24.df$ASV10)*100,
                                      sd(mix.high.T24.df$ASV2)*100,
                                      sd(mix.high.T24.df$ASV6)*100,
                                      sd(mix.high.T24.df$ASV5)*100,
                                      sd(mix.high.T24.df$ASV1)*100,
                                      sd(mix.high.T24.df$ASV3)*100,
                                      sd(mix.high.T24.df$ASV7)*100),
                               Treatment = "Mixture High",
                               Time = "24",
                               Taxon = colnames(mix.high.T24.df),
                               SamplesCombined = nrow(mix.high.T24.df))

# T48 dataframe
taxaOrder <- colnames(mix.high.T48.df)

mix.high.T48.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(mix.high.T48.df$ASV4)*100,
                                                mean(mix.high.T48.df$ASV8)*100,
                                                mean(mix.high.T48.df$ASV9)*100,
                                                mean(mix.high.T48.df$ASV10)*100,
                                                mean(mix.high.T48.df$ASV2)*100,
                                                mean(mix.high.T48.df$ASV6)*100,
                                                mean(mix.high.T48.df$ASV5)*100,
                                                mean(mix.high.T48.df$ASV1)*100,
                                                mean(mix.high.T48.df$ASV3)*100,
                                                mean(mix.high.T48.df$ASV7)*100), 
                               SD = c(sd(mix.high.T48.df$ASV4)*100,
                                      sd(mix.high.T48.df$ASV8)*100,
                                      sd(mix.high.T48.df$ASV9)*100,
                                      sd(mix.high.T48.df$ASV10)*100,
                                      sd(mix.high.T48.df$ASV2)*100,
                                      sd(mix.high.T48.df$ASV6)*100,
                                      sd(mix.high.T48.df$ASV5)*100,
                                      sd(mix.high.T48.df$ASV1)*100,
                                      sd(mix.high.T48.df$ASV3)*100,
                                      sd(mix.high.T48.df$ASV7)*100),
                               Treatment = "Mixture High",
                               Time = "48",
                               Taxon = colnames(mix.high.T48.df),
                               SamplesCombined = nrow(mix.high.T48.df))

# T96 dataframe
taxaOrder <- colnames(mix.high.T96.df)

mix.high.T96.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(mix.high.T96.df$ASV4)*100,
                                                mean(mix.high.T96.df$ASV8)*100,
                                                mean(mix.high.T96.df$ASV9)*100,
                                                mean(mix.high.T96.df$ASV10)*100,
                                                mean(mix.high.T96.df$ASV2)*100,
                                                mean(mix.high.T96.df$ASV6)*100,
                                                mean(mix.high.T96.df$ASV5)*100,
                                                mean(mix.high.T96.df$ASV1)*100,
                                                mean(mix.high.T96.df$ASV3)*100,
                                                mean(mix.high.T96.df$ASV7)*100), 
                               SD = c(sd(mix.high.T96.df$ASV4)*100,
                                      sd(mix.high.T96.df$ASV8)*100,
                                      sd(mix.high.T96.df$ASV9)*100,
                                      sd(mix.high.T96.df$ASV10)*100,
                                      sd(mix.high.T96.df$ASV2)*100,
                                      sd(mix.high.T96.df$ASV6)*100,
                                      sd(mix.high.T96.df$ASV5)*100,
                                      sd(mix.high.T96.df$ASV1)*100,
                                      sd(mix.high.T96.df$ASV3)*100,
                                      sd(mix.high.T96.df$ASV7)*100),
                               Treatment = "Mixture High",
                               Time = "96",
                               Taxon = colnames(mix.high.T96.df),
                               SamplesCombined = nrow(mix.high.T96.df))

mix.high.MRA <- bind_rows(mix.high.T0.MRA, mix.high.T12.MRA, mix.high.T24.MRA, mix.high.T48.MRA, mix.high.T96.MRA)

mix.MRA <- bind_rows(mix.low.MRA, mix.high.MRA)

readr::write_csv(mix.MRA, here::here("Data/06 - Relative Abundance - Output/03 - mixture_meanRelAbund.csv"))







## Create dataframe for mean relative abundance Ampicillin Low ====
### Subset by treatment 
ps.amp.low <- subset_samples(ps.rare.top10, (Treatment == "Ampicillin Low"))

### Subset by time 
ps.amp.low.T0 <- subset_samples(ps.amp.low, (Time == "0"))
ps.amp.low.T12 <- subset_samples(ps.amp.low, (Time == "12"))
ps.amp.low.T24 <- subset_samples(ps.amp.low, (Time == "24"))
ps.amp.low.T48 <- subset_samples(ps.amp.low, (Time == "48"))
ps.amp.low.T96 <- subset_samples(ps.amp.low, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
amp.low.T0.df <- as.data.frame(ps.amp.low.T0@otu_table)
amp.low.T12.df <- as.data.frame(ps.amp.low.T12@otu_table)
amp.low.T24.df <- as.data.frame(ps.amp.low.T24@otu_table)
amp.low.T48.df <- as.data.frame(ps.amp.low.T48@otu_table)
amp.low.T96.df <- as.data.frame(ps.amp.low.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(amp.low.T0.df)

amp.low.T0.MRA <- data.frame(ASV = taxaOrder,
                             MeanRelAbund = c(mean(amp.low.T0.df$ASV4)*100,
                                              mean(amp.low.T0.df$ASV8)*100,
                                              mean(amp.low.T0.df$ASV9)*100,
                                              mean(amp.low.T0.df$ASV10)*100,
                                              mean(amp.low.T0.df$ASV2)*100,
                                              mean(amp.low.T0.df$ASV6)*100,
                                              mean(amp.low.T0.df$ASV5)*100,
                                              mean(amp.low.T0.df$ASV1)*100,
                                              mean(amp.low.T0.df$ASV3)*100,
                                              mean(amp.low.T0.df$ASV7)*100), 
                             SD = c(sd(amp.low.T0.df$ASV4)*100,
                                    sd(amp.low.T0.df$ASV8)*100,
                                    sd(amp.low.T0.df$ASV9)*100,
                                    sd(amp.low.T0.df$ASV10)*100,
                                    sd(amp.low.T0.df$ASV2)*100,
                                    sd(amp.low.T0.df$ASV6)*100,
                                    sd(amp.low.T0.df$ASV5)*100,
                                    sd(amp.low.T0.df$ASV1)*100,
                                    sd(amp.low.T0.df$ASV3)*100,
                                    sd(amp.low.T0.df$ASV7)*100),
                             Treatment = "Ampicillin Low",
                             Time = "0",
                             Taxon = colnames(amp.low.T0.df),
                             SamplesCombined = nrow(amp.low.T0.df))

### T12 dataframe
taxaOrder <- colnames(amp.low.T12.df)

amp.low.T12.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(amp.low.T12.df$ASV4)*100,
                                               mean(amp.low.T12.df$ASV8)*100,
                                               mean(amp.low.T12.df$ASV9)*100,
                                               mean(amp.low.T12.df$ASV10)*100,
                                               mean(amp.low.T12.df$ASV2)*100,
                                               mean(amp.low.T12.df$ASV6)*100,
                                               mean(amp.low.T12.df$ASV5)*100,
                                               mean(amp.low.T12.df$ASV1)*100,
                                               mean(amp.low.T12.df$ASV3)*100,
                                               mean(amp.low.T12.df$ASV7)*100), 
                              SD = c(sd(amp.low.T12.df$ASV4)*100,
                                     sd(amp.low.T12.df$ASV8)*100,
                                     sd(amp.low.T12.df$ASV9)*100,
                                     sd(amp.low.T12.df$ASV10)*100,
                                     sd(amp.low.T12.df$ASV2)*100,
                                     sd(amp.low.T12.df$ASV6)*100,
                                     sd(amp.low.T12.df$ASV5)*100,
                                     sd(amp.low.T12.df$ASV1)*100,
                                     sd(amp.low.T12.df$ASV3)*100,
                                     sd(amp.low.T12.df$ASV7)*100),
                              Treatment = "Ampicillin Low",
                              Time = "12",
                              Taxon = colnames(amp.low.T12.df),
                              SamplesCombined = nrow(amp.low.T12.df))

# T24 dataframe
taxaOrder <- colnames(amp.low.T24.df)

amp.low.T24.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(amp.low.T24.df$ASV4)*100,
                                               mean(amp.low.T24.df$ASV8)*100,
                                               mean(amp.low.T24.df$ASV9)*100,
                                               mean(amp.low.T24.df$ASV10)*100,
                                               mean(amp.low.T24.df$ASV2)*100,
                                               mean(amp.low.T24.df$ASV6)*100,
                                               mean(amp.low.T24.df$ASV5)*100,
                                               mean(amp.low.T24.df$ASV1)*100,
                                               mean(amp.low.T24.df$ASV3)*100,
                                               mean(amp.low.T24.df$ASV7)*100), 
                              SD = c(sd(amp.low.T24.df$ASV4)*100,
                                     sd(amp.low.T24.df$ASV8)*100,
                                     sd(amp.low.T24.df$ASV9)*100,
                                     sd(amp.low.T24.df$ASV10)*100,
                                     sd(amp.low.T24.df$ASV2)*100,
                                     sd(amp.low.T24.df$ASV6)*100,
                                     sd(amp.low.T24.df$ASV5)*100,
                                     sd(amp.low.T24.df$ASV1)*100,
                                     sd(amp.low.T24.df$ASV3)*100,
                                     sd(amp.low.T24.df$ASV7)*100),
                              Treatment = "Ampicillin Low",
                              Time = "24",
                              Taxon = colnames(amp.low.T24.df),
                              SamplesCombined = nrow(amp.low.T24.df))

# T48 dataframe
taxaOrder <- colnames(amp.low.T48.df)

amp.low.T48.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(amp.low.T48.df$ASV4)*100,
                                               mean(amp.low.T48.df$ASV8)*100,
                                               mean(amp.low.T48.df$ASV9)*100,
                                               mean(amp.low.T48.df$ASV10)*100,
                                               mean(amp.low.T48.df$ASV2)*100,
                                               mean(amp.low.T48.df$ASV6)*100,
                                               mean(amp.low.T48.df$ASV5)*100,
                                               mean(amp.low.T48.df$ASV1)*100,
                                               mean(amp.low.T48.df$ASV3)*100,
                                               mean(amp.low.T48.df$ASV7)*100), 
                              SD = c(sd(amp.low.T48.df$ASV4)*100,
                                     sd(amp.low.T48.df$ASV8)*100,
                                     sd(amp.low.T48.df$ASV9)*100,
                                     sd(amp.low.T48.df$ASV10)*100,
                                     sd(amp.low.T48.df$ASV2)*100,
                                     sd(amp.low.T48.df$ASV6)*100,
                                     sd(amp.low.T48.df$ASV5)*100,
                                     sd(amp.low.T48.df$ASV1)*100,
                                     sd(amp.low.T48.df$ASV3)*100,
                                     sd(amp.low.T48.df$ASV7)*100),
                              Treatment = "Ampicillin Low",
                              Time = "48",
                              Taxon = colnames(amp.low.T48.df),
                              SamplesCombined = nrow(amp.low.T48.df))

# T96 dataframe
taxaOrder <- colnames(amp.low.T96.df)

amp.low.T96.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(amp.low.T96.df$ASV4)*100,
                                               mean(amp.low.T96.df$ASV8)*100,
                                               mean(amp.low.T96.df$ASV9)*100,
                                               mean(amp.low.T96.df$ASV10)*100,
                                               mean(amp.low.T96.df$ASV2)*100,
                                               mean(amp.low.T96.df$ASV6)*100,
                                               mean(amp.low.T96.df$ASV5)*100,
                                               mean(amp.low.T96.df$ASV1)*100,
                                               mean(amp.low.T96.df$ASV3)*100,
                                               mean(amp.low.T96.df$ASV7)*100), 
                              SD = c(sd(amp.low.T96.df$ASV4)*100,
                                     sd(amp.low.T96.df$ASV8)*100,
                                     sd(amp.low.T96.df$ASV9)*100,
                                     sd(amp.low.T96.df$ASV10)*100,
                                     sd(amp.low.T96.df$ASV2)*100,
                                     sd(amp.low.T96.df$ASV6)*100,
                                     sd(amp.low.T96.df$ASV5)*100,
                                     sd(amp.low.T96.df$ASV1)*100,
                                     sd(amp.low.T96.df$ASV3)*100,
                                     sd(amp.low.T96.df$ASV7)*100),
                              Treatment = "Ampicillin Low",
                              Time = "96",
                              Taxon = colnames(amp.low.T96.df),
                              SamplesCombined = nrow(amp.low.T96.df))

amp.low.MRA <- bind_rows(amp.low.T0.MRA, amp.low.T12.MRA, amp.low.T24.MRA, amp.low.T48.MRA, amp.low.T96.MRA)

## Create dataframe for mean relative abundance Ampicillin High ====

### Subset by treatment 
ps.amp.high <- subset_samples(ps.rare.top10, (Treatment == "Ampicillin High"))

### Subset by time 
ps.amp.high.T0 <- subset_samples(ps.amp.high, (Time == "0"))
ps.amp.high.T12 <- subset_samples(ps.amp.high, (Time == "12"))
ps.amp.high.T24 <- subset_samples(ps.amp.high, (Time == "24"))
ps.amp.high.T48 <- subset_samples(ps.amp.high, (Time == "48"))
ps.amp.high.T96 <- subset_samples(ps.amp.high, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
amp.high.T0.df <- as.data.frame(ps.amp.high.T0@otu_table)
amp.high.T12.df <- as.data.frame(ps.amp.high.T12@otu_table)
amp.high.T24.df <- as.data.frame(ps.amp.high.T24@otu_table)
amp.high.T48.df <- as.data.frame(ps.amp.high.T48@otu_table)
amp.high.T96.df <- as.data.frame(ps.amp.high.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(amp.high.T0.df)

amp.high.T0.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(amp.high.T0.df$ASV4)*100,
                                               mean(amp.high.T0.df$ASV8)*100,
                                               mean(amp.high.T0.df$ASV9)*100,
                                               mean(amp.high.T0.df$ASV10)*100,
                                               mean(amp.high.T0.df$ASV2)*100,
                                               mean(amp.high.T0.df$ASV6)*100,
                                               mean(amp.high.T0.df$ASV5)*100,
                                               mean(amp.high.T0.df$ASV1)*100,
                                               mean(amp.high.T0.df$ASV3)*100,
                                               mean(amp.high.T0.df$ASV7)*100), 
                              SD = c(sd(amp.high.T0.df$ASV4)*100,
                                     sd(amp.high.T0.df$ASV8)*100,
                                     sd(amp.high.T0.df$ASV9)*100,
                                     sd(amp.high.T0.df$ASV10)*100,
                                     sd(amp.high.T0.df$ASV2)*100,
                                     sd(amp.high.T0.df$ASV6)*100,
                                     sd(amp.high.T0.df$ASV5)*100,
                                     sd(amp.high.T0.df$ASV1)*100,
                                     sd(amp.high.T0.df$ASV3)*100,
                                     sd(amp.high.T0.df$ASV7)*100),
                              Treatment = "Ampicillin High",
                              Time = "0",
                              Taxon = colnames(amp.high.T0.df),
                              SamplesCombined = nrow(amp.high.T0.df))

### T12 dataframe
taxaOrder <- colnames(amp.high.T12.df)

amp.high.T12.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(amp.high.T12.df$ASV4)*100,
                                                mean(amp.high.T12.df$ASV8)*100,
                                                mean(amp.high.T12.df$ASV9)*100,
                                                mean(amp.high.T12.df$ASV10)*100,
                                                mean(amp.high.T12.df$ASV2)*100,
                                                mean(amp.high.T12.df$ASV6)*100,
                                                mean(amp.high.T12.df$ASV5)*100,
                                                mean(amp.high.T12.df$ASV1)*100,
                                                mean(amp.high.T12.df$ASV3)*100,
                                                mean(amp.high.T12.df$ASV7)*100), 
                               SD = c(sd(amp.high.T12.df$ASV4)*100,
                                      sd(amp.high.T12.df$ASV8)*100,
                                      sd(amp.high.T12.df$ASV9)*100,
                                      sd(amp.high.T12.df$ASV10)*100,
                                      sd(amp.high.T12.df$ASV2)*100,
                                      sd(amp.high.T12.df$ASV6)*100,
                                      sd(amp.high.T12.df$ASV5)*100,
                                      sd(amp.high.T12.df$ASV1)*100,
                                      sd(amp.high.T12.df$ASV3)*100,
                                      sd(amp.high.T12.df$ASV7)*100),
                               Treatment = "Ampicillin High",
                               Time = "12",
                               Taxon = colnames(amp.high.T12.df),
                               SamplesCombined = nrow(amp.high.T12.df))

# T24 dataframe
taxaOrder <- colnames(amp.high.T24.df)

amp.high.T24.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(amp.high.T24.df$ASV4)*100,
                                                mean(amp.high.T24.df$ASV8)*100,
                                                mean(amp.high.T24.df$ASV9)*100,
                                                mean(amp.high.T24.df$ASV10)*100,
                                                mean(amp.high.T24.df$ASV2)*100,
                                                mean(amp.high.T24.df$ASV6)*100,
                                                mean(amp.high.T24.df$ASV5)*100,
                                                mean(amp.high.T24.df$ASV1)*100,
                                                mean(amp.high.T24.df$ASV3)*100,
                                                mean(amp.high.T24.df$ASV7)*100), 
                               SD = c(sd(amp.high.T24.df$ASV4)*100,
                                      sd(amp.high.T24.df$ASV8)*100,
                                      sd(amp.high.T24.df$ASV9)*100,
                                      sd(amp.high.T24.df$ASV10)*100,
                                      sd(amp.high.T24.df$ASV2)*100,
                                      sd(amp.high.T24.df$ASV6)*100,
                                      sd(amp.high.T24.df$ASV5)*100,
                                      sd(amp.high.T24.df$ASV1)*100,
                                      sd(amp.high.T24.df$ASV3)*100,
                                      sd(amp.high.T24.df$ASV7)*100),
                               Treatment = "Ampicillin High",
                               Time = "24",
                               Taxon = colnames(amp.high.T24.df),
                               SamplesCombined = nrow(amp.high.T24.df))

# T48 dataframe
taxaOrder <- colnames(amp.high.T48.df)

amp.high.T48.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(amp.high.T48.df$ASV4)*100,
                                                mean(amp.high.T48.df$ASV8)*100,
                                                mean(amp.high.T48.df$ASV9)*100,
                                                mean(amp.high.T48.df$ASV10)*100,
                                                mean(amp.high.T48.df$ASV2)*100,
                                                mean(amp.high.T48.df$ASV6)*100,
                                                mean(amp.high.T48.df$ASV5)*100,
                                                mean(amp.high.T48.df$ASV1)*100,
                                                mean(amp.high.T48.df$ASV3)*100,
                                                mean(amp.high.T48.df$ASV7)*100), 
                               SD = c(sd(amp.high.T48.df$ASV4)*100,
                                      sd(amp.high.T48.df$ASV8)*100,
                                      sd(amp.high.T48.df$ASV9)*100,
                                      sd(amp.high.T48.df$ASV10)*100,
                                      sd(amp.high.T48.df$ASV2)*100,
                                      sd(amp.high.T48.df$ASV6)*100,
                                      sd(amp.high.T48.df$ASV5)*100,
                                      sd(amp.high.T48.df$ASV1)*100,
                                      sd(amp.high.T48.df$ASV3)*100,
                                      sd(amp.high.T48.df$ASV7)*100),
                               Treatment = "Ampicillin High",
                               Time = "48",
                               Taxon = colnames(amp.high.T48.df),
                               SamplesCombined = nrow(amp.high.T48.df))

# T96 dataframe
taxaOrder <- colnames(amp.high.T96.df)

amp.high.T96.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(amp.high.T96.df$ASV4)*100,
                                                mean(amp.high.T96.df$ASV8)*100,
                                                mean(amp.high.T96.df$ASV9)*100,
                                                mean(amp.high.T96.df$ASV10)*100,
                                                mean(amp.high.T96.df$ASV2)*100,
                                                mean(amp.high.T96.df$ASV6)*100,
                                                mean(amp.high.T96.df$ASV5)*100,
                                                mean(amp.high.T96.df$ASV1)*100,
                                                mean(amp.high.T96.df$ASV3)*100,
                                                mean(amp.high.T96.df$ASV7)*100), 
                               SD = c(sd(amp.high.T96.df$ASV4)*100,
                                      sd(amp.high.T96.df$ASV8)*100,
                                      sd(amp.high.T96.df$ASV9)*100,
                                      sd(amp.high.T96.df$ASV10)*100,
                                      sd(amp.high.T96.df$ASV2)*100,
                                      sd(amp.high.T96.df$ASV6)*100,
                                      sd(amp.high.T96.df$ASV5)*100,
                                      sd(amp.high.T96.df$ASV1)*100,
                                      sd(amp.high.T96.df$ASV3)*100,
                                      sd(amp.high.T96.df$ASV7)*100),
                               Treatment = "Ampicillin High",
                               Time = "96",
                               Taxon = colnames(amp.high.T96.df),
                               SamplesCombined = nrow(amp.high.T96.df))

amp.high.MRA <- bind_rows(amp.high.T0.MRA, amp.high.T12.MRA, amp.high.T24.MRA, amp.high.T48.MRA, amp.high.T96.MRA)

amp.MRA <- bind_rows(amp.low.MRA, amp.high.MRA)

readr::write_csv(amp.MRA, here::here("Data/06 - Relative Abundance - Output/03 - ampicillin_meanRelAbund.csv"))


## Create dataframe for mean relative abundance Streptomycin Low ====
### Subset by treatment 
ps.strep.low <- subset_samples(ps.rare.top10, (Treatment == "Streptomycin Low"))

### Subset by time 
ps.strep.low.T0 <- subset_samples(ps.strep.low, (Time == "0"))
ps.strep.low.T12 <- subset_samples(ps.strep.low, (Time == "12"))
ps.strep.low.T24 <- subset_samples(ps.strep.low, (Time == "24"))
ps.strep.low.T48 <- subset_samples(ps.strep.low, (Time == "48"))
ps.strep.low.T96 <- subset_samples(ps.strep.low, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
strep.low.T0.df <- as.data.frame(ps.strep.low.T0@otu_table)
strep.low.T12.df <- as.data.frame(ps.strep.low.T12@otu_table)
strep.low.T24.df <- as.data.frame(ps.strep.low.T24@otu_table)
strep.low.T48.df <- as.data.frame(ps.strep.low.T48@otu_table)
strep.low.T96.df <- as.data.frame(ps.strep.low.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(strep.low.T0.df)

strep.low.T0.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(strep.low.T0.df$ASV4)*100,
                                                mean(strep.low.T0.df$ASV8)*100,
                                                mean(strep.low.T0.df$ASV9)*100,
                                                mean(strep.low.T0.df$ASV10)*100,
                                                mean(strep.low.T0.df$ASV2)*100,
                                                mean(strep.low.T0.df$ASV6)*100,
                                                mean(strep.low.T0.df$ASV5)*100,
                                                mean(strep.low.T0.df$ASV1)*100,
                                                mean(strep.low.T0.df$ASV3)*100,
                                                mean(strep.low.T0.df$ASV7)*100), 
                               SD = c(sd(strep.low.T0.df$ASV4)*100,
                                      sd(strep.low.T0.df$ASV8)*100,
                                      sd(strep.low.T0.df$ASV9)*100,
                                      sd(strep.low.T0.df$ASV10)*100,
                                      sd(strep.low.T0.df$ASV2)*100,
                                      sd(strep.low.T0.df$ASV6)*100,
                                      sd(strep.low.T0.df$ASV5)*100,
                                      sd(strep.low.T0.df$ASV1)*100,
                                      sd(strep.low.T0.df$ASV3)*100,
                                      sd(strep.low.T0.df$ASV7)*100),
                               Treatment = "Streptomycin Low",
                               Time = "0",
                               Taxon = colnames(strep.low.T0.df),
                               SamplesCombined = nrow(strep.low.T0.df))

### T12 dataframe
taxaOrder <- colnames(strep.low.T12.df)

strep.low.T12.MRA <- data.frame(ASV = taxaOrder,
                                MeanRelAbund = c(mean(strep.low.T12.df$ASV4)*100,
                                                 mean(strep.low.T12.df$ASV8)*100,
                                                 mean(strep.low.T12.df$ASV9)*100,
                                                 mean(strep.low.T12.df$ASV10)*100,
                                                 mean(strep.low.T12.df$ASV2)*100,
                                                 mean(strep.low.T12.df$ASV6)*100,
                                                 mean(strep.low.T12.df$ASV5)*100,
                                                 mean(strep.low.T12.df$ASV1)*100,
                                                 mean(strep.low.T12.df$ASV3)*100,
                                                 mean(strep.low.T12.df$ASV7)*100), 
                                SD = c(sd(strep.low.T12.df$ASV4)*100,
                                       sd(strep.low.T12.df$ASV8)*100,
                                       sd(strep.low.T12.df$ASV9)*100,
                                       sd(strep.low.T12.df$ASV10)*100,
                                       sd(strep.low.T12.df$ASV2)*100,
                                       sd(strep.low.T12.df$ASV6)*100,
                                       sd(strep.low.T12.df$ASV5)*100,
                                       sd(strep.low.T12.df$ASV1)*100,
                                       sd(strep.low.T12.df$ASV3)*100,
                                       sd(strep.low.T12.df$ASV7)*100),
                                Treatment = "Streptomycin Low",
                                Time = "12",
                                Taxon = colnames(strep.low.T12.df),
                                SamplesCombined = nrow(strep.low.T12.df))

# T24 dataframe
taxaOrder <- colnames(strep.low.T24.df)

strep.low.T24.MRA <- data.frame(ASV = taxaOrder,
                                MeanRelAbund = c(mean(strep.low.T24.df$ASV4)*100,
                                                 mean(strep.low.T24.df$ASV8)*100,
                                                 mean(strep.low.T24.df$ASV9)*100,
                                                 mean(strep.low.T24.df$ASV10)*100,
                                                 mean(strep.low.T24.df$ASV2)*100,
                                                 mean(strep.low.T24.df$ASV6)*100,
                                                 mean(strep.low.T24.df$ASV5)*100,
                                                 mean(strep.low.T24.df$ASV1)*100,
                                                 mean(strep.low.T24.df$ASV3)*100,
                                                 mean(strep.low.T24.df$ASV7)*100), 
                                SD = c(sd(strep.low.T24.df$ASV4)*100,
                                       sd(strep.low.T24.df$ASV8)*100,
                                       sd(strep.low.T24.df$ASV9)*100,
                                       sd(strep.low.T24.df$ASV10)*100,
                                       sd(strep.low.T24.df$ASV2)*100,
                                       sd(strep.low.T24.df$ASV6)*100,
                                       sd(strep.low.T24.df$ASV5)*100,
                                       sd(strep.low.T24.df$ASV1)*100,
                                       sd(strep.low.T24.df$ASV3)*100,
                                       sd(strep.low.T24.df$ASV7)*100),
                                Treatment = "Streptomycin Low",
                                Time = "24",
                                Taxon = colnames(strep.low.T24.df),
                                SamplesCombined = nrow(strep.low.T24.df))

# T48 dataframe
taxaOrder <- colnames(strep.low.T48.df)

strep.low.T48.MRA <- data.frame(ASV = taxaOrder,
                                MeanRelAbund = c(mean(strep.low.T48.df$ASV4)*100,
                                                 mean(strep.low.T48.df$ASV8)*100,
                                                 mean(strep.low.T48.df$ASV9)*100,
                                                 mean(strep.low.T48.df$ASV10)*100,
                                                 mean(strep.low.T48.df$ASV2)*100,
                                                 mean(strep.low.T48.df$ASV6)*100,
                                                 mean(strep.low.T48.df$ASV5)*100,
                                                 mean(strep.low.T48.df$ASV1)*100,
                                                 mean(strep.low.T48.df$ASV3)*100,
                                                 mean(strep.low.T48.df$ASV7)*100), 
                                SD = c(sd(strep.low.T48.df$ASV4)*100,
                                       sd(strep.low.T48.df$ASV8)*100,
                                       sd(strep.low.T48.df$ASV9)*100,
                                       sd(strep.low.T48.df$ASV10)*100,
                                       sd(strep.low.T48.df$ASV2)*100,
                                       sd(strep.low.T48.df$ASV6)*100,
                                       sd(strep.low.T48.df$ASV5)*100,
                                       sd(strep.low.T48.df$ASV1)*100,
                                       sd(strep.low.T48.df$ASV3)*100,
                                       sd(strep.low.T48.df$ASV7)*100),
                                Treatment = "Streptomycin Low",
                                Time = "48",
                                Taxon = colnames(strep.low.T48.df),
                                SamplesCombined = nrow(strep.low.T48.df))

# T96 dataframe
taxaOrder <- colnames(strep.low.T96.df)

strep.low.T96.MRA <- data.frame(ASV = taxaOrder,
                                MeanRelAbund = c(mean(strep.low.T96.df$ASV4)*100,
                                                 mean(strep.low.T96.df$ASV8)*100,
                                                 mean(strep.low.T96.df$ASV9)*100,
                                                 mean(strep.low.T96.df$ASV10)*100,
                                                 mean(strep.low.T96.df$ASV2)*100,
                                                 mean(strep.low.T96.df$ASV6)*100,
                                                 mean(strep.low.T96.df$ASV5)*100,
                                                 mean(strep.low.T96.df$ASV1)*100,
                                                 mean(strep.low.T96.df$ASV3)*100,
                                                 mean(strep.low.T96.df$ASV7)*100), 
                                SD = c(sd(strep.low.T96.df$ASV4)*100,
                                       sd(strep.low.T96.df$ASV8)*100,
                                       sd(strep.low.T96.df$ASV9)*100,
                                       sd(strep.low.T96.df$ASV10)*100,
                                       sd(strep.low.T96.df$ASV2)*100,
                                       sd(strep.low.T96.df$ASV6)*100,
                                       sd(strep.low.T96.df$ASV5)*100,
                                       sd(strep.low.T96.df$ASV1)*100,
                                       sd(strep.low.T96.df$ASV3)*100,
                                       sd(strep.low.T96.df$ASV7)*100),
                                Treatment = "Streptomycin Low",
                                Time = "96",
                                Taxon = colnames(strep.low.T96.df),
                                SamplesCombined = nrow(strep.low.T96.df))

strep.low.MRA <- bind_rows(strep.low.T0.MRA, strep.low.T12.MRA, strep.low.T24.MRA, strep.low.T48.MRA, strep.low.T96.MRA)

## Create dataframe for mean relative abundance Streptomycin High ====
### Subset by treatment 
ps.strep.high <- subset_samples(ps.rare.top10, (Treatment == "Streptomycin High"))

### Subset by time 
ps.strep.high.T0 <- subset_samples(ps.strep.high, (Time == "0"))
ps.strep.high.T12 <- subset_samples(ps.strep.high, (Time == "12"))
ps.strep.high.T24 <- subset_samples(ps.strep.high, (Time == "24"))
ps.strep.high.T48 <- subset_samples(ps.strep.high, (Time == "48"))
ps.strep.high.T96 <- subset_samples(ps.strep.high, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
strep.high.T0.df <- as.data.frame(ps.strep.high.T0@otu_table)
strep.high.T12.df <- as.data.frame(ps.strep.high.T12@otu_table)
strep.high.T24.df <- as.data.frame(ps.strep.high.T24@otu_table)
strep.high.T48.df <- as.data.frame(ps.strep.high.T48@otu_table)
strep.high.T96.df <- as.data.frame(ps.strep.high.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(strep.high.T0.df)

strep.high.T0.MRA <- data.frame(ASV = taxaOrder,
                                MeanRelAbund = c(mean(strep.high.T0.df$ASV4)*100,
                                                 mean(strep.high.T0.df$ASV8)*100,
                                                 mean(strep.high.T0.df$ASV9)*100,
                                                 mean(strep.high.T0.df$ASV10)*100,
                                                 mean(strep.high.T0.df$ASV2)*100,
                                                 mean(strep.high.T0.df$ASV6)*100,
                                                 mean(strep.high.T0.df$ASV5)*100,
                                                 mean(strep.high.T0.df$ASV1)*100,
                                                 mean(strep.high.T0.df$ASV3)*100,
                                                 mean(strep.high.T0.df$ASV7)*100), 
                                SD = c(sd(strep.high.T0.df$ASV4)*100,
                                       sd(strep.high.T0.df$ASV8)*100,
                                       sd(strep.high.T0.df$ASV9)*100,
                                       sd(strep.high.T0.df$ASV10)*100,
                                       sd(strep.high.T0.df$ASV2)*100,
                                       sd(strep.high.T0.df$ASV6)*100,
                                       sd(strep.high.T0.df$ASV5)*100,
                                       sd(strep.high.T0.df$ASV1)*100,
                                       sd(strep.high.T0.df$ASV3)*100,
                                       sd(strep.high.T0.df$ASV7)*100),
                                Treatment = "Streptomycin High",
                                Time = "0",
                                Taxon = colnames(strep.high.T0.df),
                                SamplesCombined = nrow(strep.high.T0.df))

### T12 dataframe
taxaOrder <- colnames(strep.high.T12.df)

strep.high.T12.MRA <- data.frame(ASV = taxaOrder,
                                 MeanRelAbund = c(mean(strep.high.T12.df$ASV4)*100,
                                                  mean(strep.high.T12.df$ASV8)*100,
                                                  mean(strep.high.T12.df$ASV9)*100,
                                                  mean(strep.high.T12.df$ASV10)*100,
                                                  mean(strep.high.T12.df$ASV2)*100,
                                                  mean(strep.high.T12.df$ASV6)*100,
                                                  mean(strep.high.T12.df$ASV5)*100,
                                                  mean(strep.high.T12.df$ASV1)*100,
                                                  mean(strep.high.T12.df$ASV3)*100,
                                                  mean(strep.high.T12.df$ASV7)*100), 
                                 SD = c(sd(strep.high.T12.df$ASV4)*100,
                                        sd(strep.high.T12.df$ASV8)*100,
                                        sd(strep.high.T12.df$ASV9)*100,
                                        sd(strep.high.T12.df$ASV10)*100,
                                        sd(strep.high.T12.df$ASV2)*100,
                                        sd(strep.high.T12.df$ASV6)*100,
                                        sd(strep.high.T12.df$ASV5)*100,
                                        sd(strep.high.T12.df$ASV1)*100,
                                        sd(strep.high.T12.df$ASV3)*100,
                                        sd(strep.high.T12.df$ASV7)*100),
                                 Treatment = "Streptomycin High",
                                 Time = "12",
                                 Taxon = colnames(strep.high.T12.df),
                                 SamplesCombined = nrow(strep.high.T12.df))

# T24 dataframe
taxaOrder <- colnames(strep.high.T24.df)

strep.high.T24.MRA <- data.frame(ASV = taxaOrder,
                                 MeanRelAbund = c(mean(strep.high.T24.df$ASV4)*100,
                                                  mean(strep.high.T24.df$ASV8)*100,
                                                  mean(strep.high.T24.df$ASV9)*100,
                                                  mean(strep.high.T24.df$ASV10)*100,
                                                  mean(strep.high.T24.df$ASV2)*100,
                                                  mean(strep.high.T24.df$ASV6)*100,
                                                  mean(strep.high.T24.df$ASV5)*100,
                                                  mean(strep.high.T24.df$ASV1)*100,
                                                  mean(strep.high.T24.df$ASV3)*100,
                                                  mean(strep.high.T24.df$ASV7)*100), 
                                 SD = c(sd(strep.high.T24.df$ASV4)*100,
                                        sd(strep.high.T24.df$ASV8)*100,
                                        sd(strep.high.T24.df$ASV9)*100,
                                        sd(strep.high.T24.df$ASV10)*100,
                                        sd(strep.high.T24.df$ASV2)*100,
                                        sd(strep.high.T24.df$ASV6)*100,
                                        sd(strep.high.T24.df$ASV5)*100,
                                        sd(strep.high.T24.df$ASV1)*100,
                                        sd(strep.high.T24.df$ASV3)*100,
                                        sd(strep.high.T24.df$ASV7)*100),
                                 Treatment = "Streptomycin High",
                                 Time = "24",
                                 Taxon = colnames(strep.high.T24.df),
                                 SamplesCombined = nrow(strep.high.T24.df))

# T48 dataframe
taxaOrder <- colnames(strep.high.T48.df)

strep.high.T48.MRA <- data.frame(ASV = taxaOrder,
                                 MeanRelAbund = c(mean(strep.high.T48.df$ASV4)*100,
                                                  mean(strep.high.T48.df$ASV8)*100,
                                                  mean(strep.high.T48.df$ASV9)*100,
                                                  mean(strep.high.T48.df$ASV10)*100,
                                                  mean(strep.high.T48.df$ASV2)*100,
                                                  mean(strep.high.T48.df$ASV6)*100,
                                                  mean(strep.high.T48.df$ASV5)*100,
                                                  mean(strep.high.T48.df$ASV1)*100,
                                                  mean(strep.high.T48.df$ASV3)*100,
                                                  mean(strep.high.T48.df$ASV7)*100), 
                                 SD = c(sd(strep.high.T48.df$ASV4)*100,
                                        sd(strep.high.T48.df$ASV8)*100,
                                        sd(strep.high.T48.df$ASV9)*100,
                                        sd(strep.high.T48.df$ASV10)*100,
                                        sd(strep.high.T48.df$ASV2)*100,
                                        sd(strep.high.T48.df$ASV6)*100,
                                        sd(strep.high.T48.df$ASV5)*100,
                                        sd(strep.high.T48.df$ASV1)*100,
                                        sd(strep.high.T48.df$ASV3)*100,
                                        sd(strep.high.T48.df$ASV7)*100),
                                 Treatment = "Streptomycin High",
                                 Time = "48",
                                 Taxon = colnames(strep.high.T48.df),
                                 SamplesCombined = nrow(strep.high.T48.df))

# T96 dataframe
taxaOrder <- colnames(strep.high.T96.df)

strep.high.T96.MRA <- data.frame(ASV = taxaOrder,
                                 MeanRelAbund = c(mean(strep.high.T96.df$ASV4)*100,
                                                  mean(strep.high.T96.df$ASV8)*100,
                                                  mean(strep.high.T96.df$ASV9)*100,
                                                  mean(strep.high.T96.df$ASV10)*100,
                                                  mean(strep.high.T96.df$ASV2)*100,
                                                  mean(strep.high.T96.df$ASV6)*100,
                                                  mean(strep.high.T96.df$ASV5)*100,
                                                  mean(strep.high.T96.df$ASV1)*100,
                                                  mean(strep.high.T96.df$ASV3)*100,
                                                  mean(strep.high.T96.df$ASV7)*100), 
                                 SD = c(sd(strep.high.T96.df$ASV4)*100,
                                        sd(strep.high.T96.df$ASV8)*100,
                                        sd(strep.high.T96.df$ASV9)*100,
                                        sd(strep.high.T96.df$ASV10)*100,
                                        sd(strep.high.T96.df$ASV2)*100,
                                        sd(strep.high.T96.df$ASV6)*100,
                                        sd(strep.high.T96.df$ASV5)*100,
                                        sd(strep.high.T96.df$ASV1)*100,
                                        sd(strep.high.T96.df$ASV3)*100,
                                        sd(strep.high.T96.df$ASV7)*100),
                                 Treatment = "Streptomycin High",
                                 Time = "96",
                                 Taxon = colnames(strep.high.T96.df),
                                 SamplesCombined = nrow(strep.high.T96.df))

strep.high.MRA <- bind_rows(strep.high.T0.MRA, strep.high.T12.MRA, strep.high.T24.MRA, strep.high.T48.MRA, strep.high.T96.MRA)

strep.MRA <- bind_rows(strep.low.MRA, strep.high.MRA)

readr::write_csv(strep.MRA, here::here("Data/06 - Relative Abundance - Output/04 - streptomycin_meanRelAbund.csv"))


## Create dataframe for mean relative abundance Ciprofloxacin Low ====
### Subset by treatment 
ps.cip.low <- subset_samples(ps.rare.top10, (Treatment == "Ciprofloxacin Low"))

### Subset by time 
ps.cip.low.T0 <- subset_samples(ps.cip.low, (Time == "0"))
ps.cip.low.T12 <- subset_samples(ps.cip.low, (Time == "12"))
ps.cip.low.T24 <- subset_samples(ps.cip.low, (Time == "24"))
ps.cip.low.T48 <- subset_samples(ps.cip.low, (Time == "48"))
ps.cip.low.T96 <- subset_samples(ps.cip.low, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
cip.low.T0.df <- as.data.frame(ps.cip.low.T0@otu_table)
cip.low.T12.df <- as.data.frame(ps.cip.low.T12@otu_table)
cip.low.T24.df <- as.data.frame(ps.cip.low.T24@otu_table)
cip.low.T48.df <- as.data.frame(ps.cip.low.T48@otu_table)
cip.low.T96.df <- as.data.frame(ps.cip.low.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(cip.low.T0.df)

cip.low.T0.MRA <- data.frame(ASV = taxaOrder,
                             MeanRelAbund = c(mean(cip.low.T0.df$ASV4)*100,
                                              mean(cip.low.T0.df$ASV8)*100,
                                              mean(cip.low.T0.df$ASV9)*100,
                                              mean(cip.low.T0.df$ASV10)*100,
                                              mean(cip.low.T0.df$ASV2)*100,
                                              mean(cip.low.T0.df$ASV6)*100,
                                              mean(cip.low.T0.df$ASV5)*100,
                                              mean(cip.low.T0.df$ASV1)*100,
                                              mean(cip.low.T0.df$ASV3)*100,
                                              mean(cip.low.T0.df$ASV7)*100), 
                             SD = c(sd(cip.low.T0.df$ASV4)*100,
                                    sd(cip.low.T0.df$ASV8)*100,
                                    sd(cip.low.T0.df$ASV9)*100,
                                    sd(cip.low.T0.df$ASV10)*100,
                                    sd(cip.low.T0.df$ASV2)*100,
                                    sd(cip.low.T0.df$ASV6)*100,
                                    sd(cip.low.T0.df$ASV5)*100,
                                    sd(cip.low.T0.df$ASV1)*100,
                                    sd(cip.low.T0.df$ASV3)*100,
                                    sd(cip.low.T0.df$ASV7)*100),
                             Treatment = "Ciprofloxacin Low",
                             Time = "0",
                             Taxon = colnames(cip.low.T0.df),
                             SamplesCombined = nrow(cip.low.T0.df))

### T12 dataframe
taxaOrder <- colnames(cip.low.T12.df)

cip.low.T12.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(cip.low.T12.df$ASV4)*100,
                                               mean(cip.low.T12.df$ASV8)*100,
                                               mean(cip.low.T12.df$ASV9)*100,
                                               mean(cip.low.T12.df$ASV10)*100,
                                               mean(cip.low.T12.df$ASV2)*100,
                                               mean(cip.low.T12.df$ASV6)*100,
                                               mean(cip.low.T12.df$ASV5)*100,
                                               mean(cip.low.T12.df$ASV1)*100,
                                               mean(cip.low.T12.df$ASV3)*100,
                                               mean(cip.low.T12.df$ASV7)*100), 
                              SD = c(sd(cip.low.T12.df$ASV4)*100,
                                     sd(cip.low.T12.df$ASV8)*100,
                                     sd(cip.low.T12.df$ASV9)*100,
                                     sd(cip.low.T12.df$ASV10)*100,
                                     sd(cip.low.T12.df$ASV2)*100,
                                     sd(cip.low.T12.df$ASV6)*100,
                                     sd(cip.low.T12.df$ASV5)*100,
                                     sd(cip.low.T12.df$ASV1)*100,
                                     sd(cip.low.T12.df$ASV3)*100,
                                     sd(cip.low.T12.df$ASV7)*100),
                              Treatment = "Ciprofloxacin Low",
                              Time = "12",
                              Taxon = colnames(cip.low.T12.df),
                              SamplesCombined = nrow(cip.low.T12.df))

# T24 dataframe
taxaOrder <- colnames(cip.low.T24.df)

cip.low.T24.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(cip.low.T24.df$ASV4)*100,
                                               mean(cip.low.T24.df$ASV8)*100,
                                               mean(cip.low.T24.df$ASV9)*100,
                                               mean(cip.low.T24.df$ASV10)*100,
                                               mean(cip.low.T24.df$ASV2)*100,
                                               mean(cip.low.T24.df$ASV6)*100,
                                               mean(cip.low.T24.df$ASV5)*100,
                                               mean(cip.low.T24.df$ASV1)*100,
                                               mean(cip.low.T24.df$ASV3)*100,
                                               mean(cip.low.T24.df$ASV7)*100), 
                              SD = c(sd(cip.low.T24.df$ASV4)*100,
                                     sd(cip.low.T24.df$ASV8)*100,
                                     sd(cip.low.T24.df$ASV9)*100,
                                     sd(cip.low.T24.df$ASV10)*100,
                                     sd(cip.low.T24.df$ASV2)*100,
                                     sd(cip.low.T24.df$ASV6)*100,
                                     sd(cip.low.T24.df$ASV5)*100,
                                     sd(cip.low.T24.df$ASV1)*100,
                                     sd(cip.low.T24.df$ASV3)*100,
                                     sd(cip.low.T24.df$ASV7)*100),
                              Treatment = "Ciprofloxacin Low",
                              Time = "24",
                              Taxon = colnames(cip.low.T24.df),
                              SamplesCombined = nrow(cip.low.T24.df))

# T48 dataframe
taxaOrder <- colnames(cip.low.T48.df)

cip.low.T48.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(cip.low.T48.df$ASV4)*100,
                                               mean(cip.low.T48.df$ASV8)*100,
                                               mean(cip.low.T48.df$ASV9)*100,
                                               mean(cip.low.T48.df$ASV10)*100,
                                               mean(cip.low.T48.df$ASV2)*100,
                                               mean(cip.low.T48.df$ASV6)*100,
                                               mean(cip.low.T48.df$ASV5)*100,
                                               mean(cip.low.T48.df$ASV1)*100,
                                               mean(cip.low.T48.df$ASV3)*100,
                                               mean(cip.low.T48.df$ASV7)*100), 
                              SD = c(sd(cip.low.T48.df$ASV4)*100,
                                     sd(cip.low.T48.df$ASV8)*100,
                                     sd(cip.low.T48.df$ASV9)*100,
                                     sd(cip.low.T48.df$ASV10)*100,
                                     sd(cip.low.T48.df$ASV2)*100,
                                     sd(cip.low.T48.df$ASV6)*100,
                                     sd(cip.low.T48.df$ASV5)*100,
                                     sd(cip.low.T48.df$ASV1)*100,
                                     sd(cip.low.T48.df$ASV3)*100,
                                     sd(cip.low.T48.df$ASV7)*100),
                              Treatment = "Ciprofloxacin Low",
                              Time = "48",
                              Taxon = colnames(cip.low.T48.df),
                              SamplesCombined = nrow(cip.low.T48.df))

# T96 dataframe
taxaOrder <- colnames(cip.low.T96.df)

cip.low.T96.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(cip.low.T96.df$ASV4)*100,
                                               mean(cip.low.T96.df$ASV8)*100,
                                               mean(cip.low.T96.df$ASV9)*100,
                                               mean(cip.low.T96.df$ASV10)*100,
                                               mean(cip.low.T96.df$ASV2)*100,
                                               mean(cip.low.T96.df$ASV6)*100,
                                               mean(cip.low.T96.df$ASV5)*100,
                                               mean(cip.low.T96.df$ASV1)*100,
                                               mean(cip.low.T96.df$ASV3)*100,
                                               mean(cip.low.T96.df$ASV7)*100), 
                              SD = c(sd(cip.low.T96.df$ASV4)*100,
                                     sd(cip.low.T96.df$ASV8)*100,
                                     sd(cip.low.T96.df$ASV9)*100,
                                     sd(cip.low.T96.df$ASV10)*100,
                                     sd(cip.low.T96.df$ASV2)*100,
                                     sd(cip.low.T96.df$ASV6)*100,
                                     sd(cip.low.T96.df$ASV5)*100,
                                     sd(cip.low.T96.df$ASV1)*100,
                                     sd(cip.low.T96.df$ASV3)*100,
                                     sd(cip.low.T96.df$ASV7)*100),
                              Treatment = "Ciprofloxacin Low",
                              Time = "96",
                              Taxon = colnames(cip.low.T96.df),
                              SamplesCombined = nrow(cip.low.T96.df))

cip.low.MRA <- bind_rows(cip.low.T0.MRA, cip.low.T12.MRA, cip.low.T24.MRA, cip.low.T48.MRA, cip.low.T96.MRA)

## Create dataframe for mean relative abundance Ciprofloxacin High====

### Subset by treatment 
ps.cip.high <- subset_samples(ps.rare.top10, (Treatment == "Ciprofloxacin High"))

### Subset by time 
ps.cip.high.T0 <- subset_samples(ps.cip.high, (Time == "0"))
ps.cip.high.T12 <- subset_samples(ps.cip.high, (Time == "12"))
ps.cip.high.T24 <- subset_samples(ps.cip.high, (Time == "24"))
ps.cip.high.T48 <- subset_samples(ps.cip.high, (Time == "48"))
ps.cip.high.T96 <- subset_samples(ps.cip.high, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
cip.high.T0.df <- as.data.frame(ps.cip.high.T0@otu_table)
cip.high.T12.df <- as.data.frame(ps.cip.high.T12@otu_table)
cip.high.T24.df <- as.data.frame(ps.cip.high.T24@otu_table)
cip.high.T48.df <- as.data.frame(ps.cip.high.T48@otu_table)
cip.high.T96.df <- as.data.frame(ps.cip.high.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(cip.high.T0.df)

cip.high.T0.MRA <- data.frame(ASV = taxaOrder,
                              MeanRelAbund = c(mean(cip.high.T0.df$ASV4)*100,
                                               mean(cip.high.T0.df$ASV8)*100,
                                               mean(cip.high.T0.df$ASV9)*100,
                                               mean(cip.high.T0.df$ASV10)*100,
                                               mean(cip.high.T0.df$ASV2)*100,
                                               mean(cip.high.T0.df$ASV6)*100,
                                               mean(cip.high.T0.df$ASV5)*100,
                                               mean(cip.high.T0.df$ASV1)*100,
                                               mean(cip.high.T0.df$ASV3)*100,
                                               mean(cip.high.T0.df$ASV7)*100), 
                              SD = c(sd(cip.high.T0.df$ASV4)*100,
                                     sd(cip.high.T0.df$ASV8)*100,
                                     sd(cip.high.T0.df$ASV9)*100,
                                     sd(cip.high.T0.df$ASV10)*100,
                                     sd(cip.high.T0.df$ASV2)*100,
                                     sd(cip.high.T0.df$ASV6)*100,
                                     sd(cip.high.T0.df$ASV5)*100,
                                     sd(cip.high.T0.df$ASV1)*100,
                                     sd(cip.high.T0.df$ASV3)*100,
                                     sd(cip.high.T0.df$ASV7)*100),
                              Treatment = "Ciprofloxacin High",
                              Time = "0",
                              Taxon = colnames(cip.high.T0.df),
                              SamplesCombined = nrow(cip.high.T0.df))

### T12 dataframe
taxaOrder <- colnames(cip.high.T12.df)

cip.high.T12.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(cip.high.T12.df$ASV4)*100,
                                                mean(cip.high.T12.df$ASV8)*100,
                                                mean(cip.high.T12.df$ASV9)*100,
                                                mean(cip.high.T12.df$ASV10)*100,
                                                mean(cip.high.T12.df$ASV2)*100,
                                                mean(cip.high.T12.df$ASV6)*100,
                                                mean(cip.high.T12.df$ASV5)*100,
                                                mean(cip.high.T12.df$ASV1)*100,
                                                mean(cip.high.T12.df$ASV3)*100,
                                                mean(cip.high.T12.df$ASV7)*100), 
                               SD = c(sd(cip.high.T12.df$ASV4)*100,
                                      sd(cip.high.T12.df$ASV8)*100,
                                      sd(cip.high.T12.df$ASV9)*100,
                                      sd(cip.high.T12.df$ASV10)*100,
                                      sd(cip.high.T12.df$ASV2)*100,
                                      sd(cip.high.T12.df$ASV6)*100,
                                      sd(cip.high.T12.df$ASV5)*100,
                                      sd(cip.high.T12.df$ASV1)*100,
                                      sd(cip.high.T12.df$ASV3)*100,
                                      sd(cip.high.T12.df$ASV7)*100),
                               Treatment = "Ciprofloxacin High",
                               Time = "12",
                               Taxon = colnames(cip.high.T12.df),
                               SamplesCombined = nrow(cip.high.T12.df))

# T24 dataframe
taxaOrder <- colnames(cip.high.T24.df)

cip.high.T24.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(cip.high.T24.df$ASV4)*100,
                                                mean(cip.high.T24.df$ASV8)*100,
                                                mean(cip.high.T24.df$ASV9)*100,
                                                mean(cip.high.T24.df$ASV10)*100,
                                                mean(cip.high.T24.df$ASV2)*100,
                                                mean(cip.high.T24.df$ASV6)*100,
                                                mean(cip.high.T24.df$ASV5)*100,
                                                mean(cip.high.T24.df$ASV1)*100,
                                                mean(cip.high.T24.df$ASV3)*100,
                                                mean(cip.high.T24.df$ASV7)*100), 
                               SD = c(sd(cip.high.T24.df$ASV4)*100,
                                      sd(cip.high.T24.df$ASV8)*100,
                                      sd(cip.high.T24.df$ASV9)*100,
                                      sd(cip.high.T24.df$ASV10)*100,
                                      sd(cip.high.T24.df$ASV2)*100,
                                      sd(cip.high.T24.df$ASV6)*100,
                                      sd(cip.high.T24.df$ASV5)*100,
                                      sd(cip.high.T24.df$ASV1)*100,
                                      sd(cip.high.T24.df$ASV3)*100,
                                      sd(cip.high.T24.df$ASV7)*100),
                               Treatment = "Ciprofloxacin High",
                               Time = "24",
                               Taxon = colnames(cip.high.T24.df),
                               SamplesCombined = nrow(cip.high.T24.df))

# T48 dataframe
taxaOrder <- colnames(cip.high.T48.df)

cip.high.T48.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(cip.high.T48.df$ASV4)*100,
                                                mean(cip.high.T48.df$ASV8)*100,
                                                mean(cip.high.T48.df$ASV9)*100,
                                                mean(cip.high.T48.df$ASV10)*100,
                                                mean(cip.high.T48.df$ASV2)*100,
                                                mean(cip.high.T48.df$ASV6)*100,
                                                mean(cip.high.T48.df$ASV5)*100,
                                                mean(cip.high.T48.df$ASV1)*100,
                                                mean(cip.high.T48.df$ASV3)*100,
                                                mean(cip.high.T48.df$ASV7)*100), 
                               SD = c(sd(cip.high.T48.df$ASV4)*100,
                                      sd(cip.high.T48.df$ASV8)*100,
                                      sd(cip.high.T48.df$ASV9)*100,
                                      sd(cip.high.T48.df$ASV10)*100,
                                      sd(cip.high.T48.df$ASV2)*100,
                                      sd(cip.high.T48.df$ASV6)*100,
                                      sd(cip.high.T48.df$ASV5)*100,
                                      sd(cip.high.T48.df$ASV1)*100,
                                      sd(cip.high.T48.df$ASV3)*100,
                                      sd(cip.high.T48.df$ASV7)*100),
                               Treatment = "Ciprofloxacin High",
                               Time = "48",
                               Taxon = colnames(cip.high.T48.df),
                               SamplesCombined = nrow(cip.high.T48.df))

# T96 dataframe
taxaOrder <- colnames(cip.high.T96.df)

cip.high.T96.MRA <- data.frame(ASV = taxaOrder,
                               MeanRelAbund = c(mean(cip.high.T96.df$ASV4)*100,
                                                mean(cip.high.T96.df$ASV8)*100,
                                                mean(cip.high.T96.df$ASV9)*100,
                                                mean(cip.high.T96.df$ASV10)*100,
                                                mean(cip.high.T96.df$ASV2)*100,
                                                mean(cip.high.T96.df$ASV6)*100,
                                                mean(cip.high.T96.df$ASV5)*100,
                                                mean(cip.high.T96.df$ASV1)*100,
                                                mean(cip.high.T96.df$ASV3)*100,
                                                mean(cip.high.T96.df$ASV7)*100), 
                               SD = c(sd(cip.high.T96.df$ASV4)*100,
                                      sd(cip.high.T96.df$ASV8)*100,
                                      sd(cip.high.T96.df$ASV9)*100,
                                      sd(cip.high.T96.df$ASV10)*100,
                                      sd(cip.high.T96.df$ASV2)*100,
                                      sd(cip.high.T96.df$ASV6)*100,
                                      sd(cip.high.T96.df$ASV5)*100,
                                      sd(cip.high.T96.df$ASV1)*100,
                                      sd(cip.high.T96.df$ASV3)*100,
                                      sd(cip.high.T96.df$ASV7)*100),
                               Treatment = "Ciprofloxacin High",
                               Time = "96",
                               Taxon = colnames(cip.high.T96.df),
                               SamplesCombined = nrow(cip.high.T96.df))

cip.high.MRA <- bind_rows(cip.high.T0.MRA, cip.high.T12.MRA, cip.high.T24.MRA, cip.high.T48.MRA, cip.high.T96.MRA)

cip.MRA <- bind_rows(cip.low.MRA, cip.high.MRA)

readr::write_csv(cip.MRA, here::here("Data/06 - Relative Abundance - Output/05 - ciprofloxacin_meanRelAbund.csv"))


## Create dataframe for mean relative abundance Blank ====
### Subset by treatment 
ps.blank <- subset_samples(ps.rare.top10, (Treatment == "Blank"))

### Subset by time 
ps.blank.T0 <- subset_samples(ps.blank, (Time == "0"))
ps.blank.T12 <- subset_samples(ps.blank, (Time == "12"))
ps.blank.T24 <- subset_samples(ps.blank, (Time == "24"))
ps.blank.T48 <- subset_samples(ps.blank, (Time == "48"))
ps.blank.T96 <- subset_samples(ps.blank, (Time == "96"))

### Save phyloseq object subsets otu_table at dataframe
blank.T0.df <- as.data.frame(ps.blank.T0@otu_table)
blank.T12.df <- as.data.frame(ps.blank.T12@otu_table)
blank.T24.df <- as.data.frame(ps.blank.T24@otu_table)
blank.T48.df <- as.data.frame(ps.blank.T48@otu_table)
blank.T96.df <- as.data.frame(ps.blank.T96@otu_table)

### T0 dataframe
taxaOrder <- colnames(blank.T0.df)

blank.T0.MRA <- data.frame(ASV = taxaOrder,
                           MeanRelAbund = c(mean(blank.T0.df$ASV4)*100,
                                            mean(blank.T0.df$ASV8)*100,
                                            mean(blank.T0.df$ASV9)*100,
                                            mean(blank.T0.df$ASV10)*100,
                                            mean(blank.T0.df$ASV2)*100,
                                            mean(blank.T0.df$ASV6)*100,
                                            mean(blank.T0.df$ASV5)*100,
                                            mean(blank.T0.df$ASV1)*100,
                                            mean(blank.T0.df$ASV3)*100,
                                            mean(blank.T0.df$ASV7)*100), 
                           SD = c(sd(blank.T0.df$ASV4)*100,
                                  sd(blank.T0.df$ASV8)*100,
                                  sd(blank.T0.df$ASV9)*100,
                                  sd(blank.T0.df$ASV10)*100,
                                  sd(blank.T0.df$ASV2)*100,
                                  sd(blank.T0.df$ASV6)*100,
                                  sd(blank.T0.df$ASV5)*100,
                                  sd(blank.T0.df$ASV1)*100,
                                  sd(blank.T0.df$ASV3)*100,
                                  sd(blank.T0.df$ASV7)*100),
                           Treatment = "Blank",
                           Time = "0",
                           Taxon = colnames(blank.T0.df),
                           SamplesCombined = nrow(blank.T0.df))

### T12 dataframe
taxaOrder <- colnames(blank.T12.df)

blank.T12.MRA <- data.frame(ASV = taxaOrder,
                            MeanRelAbund = c(mean(blank.T12.df$ASV4)*100,
                                             mean(blank.T12.df$ASV8)*100,
                                             mean(blank.T12.df$ASV9)*100,
                                             mean(blank.T12.df$ASV10)*100,
                                             mean(blank.T12.df$ASV2)*100,
                                             mean(blank.T12.df$ASV6)*100,
                                             mean(blank.T12.df$ASV5)*100,
                                             mean(blank.T12.df$ASV1)*100,
                                             mean(blank.T12.df$ASV3)*100,
                                             mean(blank.T12.df$ASV7)*100), 
                            SD = c(sd(blank.T12.df$ASV4)*100,
                                   sd(blank.T12.df$ASV8)*100,
                                   sd(blank.T12.df$ASV9)*100,
                                   sd(blank.T12.df$ASV10)*100,
                                   sd(blank.T12.df$ASV2)*100,
                                   sd(blank.T12.df$ASV6)*100,
                                   sd(blank.T12.df$ASV5)*100,
                                   sd(blank.T12.df$ASV1)*100,
                                   sd(blank.T12.df$ASV3)*100,
                                   sd(blank.T12.df$ASV7)*100),
                            Treatment = "Blank",
                            Time = "12",
                            Taxon = colnames(blank.T12.df),
                            SamplesCombined = nrow(blank.T12.df))

# T24 dataframe
taxaOrder <- colnames(blank.T24.df)

blank.T24.MRA <- data.frame(ASV = taxaOrder,
                            MeanRelAbund = c(mean(blank.T24.df$ASV4)*100,
                                             mean(blank.T24.df$ASV8)*100,
                                             mean(blank.T24.df$ASV9)*100,
                                             mean(blank.T24.df$ASV10)*100,
                                             mean(blank.T24.df$ASV2)*100,
                                             mean(blank.T24.df$ASV6)*100,
                                             mean(blank.T24.df$ASV5)*100,
                                             mean(blank.T24.df$ASV1)*100,
                                             mean(blank.T24.df$ASV3)*100,
                                             mean(blank.T24.df$ASV7)*100), 
                            SD = c(sd(blank.T24.df$ASV4)*100,
                                   sd(blank.T24.df$ASV8)*100,
                                   sd(blank.T24.df$ASV9)*100,
                                   sd(blank.T24.df$ASV10)*100,
                                   sd(blank.T24.df$ASV2)*100,
                                   sd(blank.T24.df$ASV6)*100,
                                   sd(blank.T24.df$ASV5)*100,
                                   sd(blank.T24.df$ASV1)*100,
                                   sd(blank.T24.df$ASV3)*100,
                                   sd(blank.T24.df$ASV7)*100),
                            Treatment = "Blank",
                            Time = "24",
                            Taxon = colnames(blank.T24.df),
                            SamplesCombined = nrow(blank.T24.df))

# T48 dataframe
taxaOrder <- colnames(blank.T48.df)

blank.T48.MRA <- data.frame(ASV = taxaOrder,
                            MeanRelAbund = c(mean(blank.T48.df$ASV4)*100,
                                             mean(blank.T48.df$ASV8)*100,
                                             mean(blank.T48.df$ASV9)*100,
                                             mean(blank.T48.df$ASV10)*100,
                                             mean(blank.T48.df$ASV2)*100,
                                             mean(blank.T48.df$ASV6)*100,
                                             mean(blank.T48.df$ASV5)*100,
                                             mean(blank.T48.df$ASV1)*100,
                                             mean(blank.T48.df$ASV3)*100,
                                             mean(blank.T48.df$ASV7)*100), 
                            SD = c(sd(blank.T48.df$ASV4)*100,
                                   sd(blank.T48.df$ASV8)*100,
                                   sd(blank.T48.df$ASV9)*100,
                                   sd(blank.T48.df$ASV10)*100,
                                   sd(blank.T48.df$ASV2)*100,
                                   sd(blank.T48.df$ASV6)*100,
                                   sd(blank.T48.df$ASV5)*100,
                                   sd(blank.T48.df$ASV1)*100,
                                   sd(blank.T48.df$ASV3)*100,
                                   sd(blank.T48.df$ASV7)*100),
                            Treatment = "Blank",
                            Time = "48",
                            Taxon = colnames(blank.T48.df),
                            SamplesCombined = nrow(blank.T48.df))

# T96 dataframe
taxaOrder <- colnames(blank.T96.df)

blank.T96.MRA <- data.frame(ASV = taxaOrder,
                            MeanRelAbund = c(mean(blank.T96.df$ASV4)*100,
                                             mean(blank.T96.df$ASV8)*100,
                                             mean(blank.T96.df$ASV9)*100,
                                             mean(blank.T96.df$ASV10)*100,
                                             mean(blank.T96.df$ASV2)*100,
                                             mean(blank.T96.df$ASV6)*100,
                                             mean(blank.T96.df$ASV5)*100,
                                             mean(blank.T96.df$ASV1)*100,
                                             mean(blank.T96.df$ASV3)*100,
                                             mean(blank.T96.df$ASV7)*100), 
                            SD = c(sd(blank.T96.df$ASV4)*100,
                                   sd(blank.T96.df$ASV8)*100,
                                   sd(blank.T96.df$ASV9)*100,
                                   sd(blank.T96.df$ASV10)*100,
                                   sd(blank.T96.df$ASV2)*100,
                                   sd(blank.T96.df$ASV6)*100,
                                   sd(blank.T96.df$ASV5)*100,
                                   sd(blank.T96.df$ASV1)*100,
                                   sd(blank.T96.df$ASV3)*100,
                                   sd(blank.T96.df$ASV7)*100),
                            Treatment = "Blank",
                            Time = "96",
                            Taxon = colnames(blank.T96.df),
                            SamplesCombined = nrow(blank.T96.df))

blank.MRA <- bind_rows(blank.T0.MRA, blank.T12.MRA, blank.T24.MRA, blank.T48.MRA, blank.T96.MRA)

readr::write_csv(blank.MRA, here::here("Data/06 - Relative Abundance - Output/06 - blank_meanRelAbund.csv"))


## Combine all dataframes and change taxon names ====
MRA_all <- bind_rows(blank.MRA, mix.MRA, amp.MRA, strep.MRA, cip.MRA)

readr::write_csv(MRA.all, here::here("Data/06 - Relative Abundance - Output/all_meanRelAbund.csv"))


# Change taxon from ASV# to actual taxon name 
MRA_all$Taxon[MRA_all$Taxon == "ASV1"] <- "Campylobacterales Order (ASV1)"
MRA_all$Taxon[MRA_all$Taxon == "ASV2"] <- "Alteromonas (ASV2)"
MRA_all$Taxon[MRA_all$Taxon == "ASV3"] <- "Helicobacteraceae Family (ASV3)"
MRA_all$Taxon[MRA_all$Taxon == "ASV4"] <- "Phaeobacter (ASV4)"
MRA_all$Taxon[MRA_all$Taxon == "ASV5"] <- "P3OB-42 (ASV5)"
MRA_all$Taxon[MRA_all$Taxon == "ASV6"] <- "Aquarickettsia (ASV6)"
MRA_all$Taxon[MRA_all$Taxon == "ASV7"] <- "Staphylococcus (ASV7)"
MRA_all$Taxon[MRA_all$Taxon == "ASV8"] <- "Caulobacteraceae Family (ASV8)"
MRA_all$Taxon[MRA_all$Taxon == "ASV9"] <- "[Caedibacter] taeniospiralis group (ASV9)"
MRA_all$Taxon[MRA_all$Taxon == "ASV10"] <- "[Caedibacter] taeniospiralis group (ASV10)"

## Plot mean relative abundance ====
MRA_all$Treatment <- factor(MRA_all$Treatment, c("Blank", "Mixture Low", "Mixture High", 
                                                 "Ampicillin Low", "Streptomycin Low","Ciprofloxacin Low",
                                                 "Ampicillin High", "Streptomycin High", "Ciprofloxacin High"), ordered = TRUE)

## Main Figure 1
plot_meanRelAbund <- ggplot(data = MRA_all, aes(x = Time, y = MeanRelAbund, fill = Taxon)) +
  scale_fill_manual(values = c("#644194", "#CAB2D6",
                               "#4A77AF", "#B3CDE1",
                               "#5D9D40", "#BDDD94",
                               "#E31A1C", "#FB9A99",
                               "#FF7F00", "#FDBF6F")) + geom_bar(stat = "identity") + 
  facet_wrap(~Treatment)

plot_meanRelAbund <- plot_meanRelAbund + labs(title = "Mean Relative Abundance") + 
  xlab("Time") + ylab("Mean Relative Abundance (%)") + guides(fill = guide_legend(title = "Genus")) + 
  theme_bw() + theme(axis.text = element_text(face = "bold", size = 11.5), axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold"))


ggplot2::ggsave(here::here("Data/06 - Relative Abundance - Output/07 - plot_meanrelabund_all.png"), plot_meanRelAbund,
                height = 450, width = 550, units = "mm",
                scale = 0.5, dpi = 1000)