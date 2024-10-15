# Author: Sunni Patton
# Last edited: 1/30/24
# Title: Running decontam

# Code largely adapted from https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

## Load libraries ====
library(decontam)
library(ggplot2)
library(phyloseq)

## Set seed ===
set.seed(123)

## Plot sample library size ====
df <- as.data.frame(sample_data(ps.All)) 
df$LibrarySize <- sample_sums(ps.All)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot_decontam <- ggplot(data=df, aes(x = Index, y = LibrarySize, color = Sample_control)) + geom_point()

## Save plot ====
ggplot2::ggsave(here::here("Data/03 - Decontam - Output/plot_decontam.png"), plot_decontam,
               height = 200, width = 350, units = "mm",
               scale = 0.5, dpi = 1000)

## Run decontam ====
sample_data(ps.All)$is.neg <- sample_data(ps.All)$Sample_control == "Negative Control"

contamdf_combined <- isContaminant(ps.All, method="combined", neg="is.neg", conc="quant_reading", threshold=0.5)
table(contamdf_combined$contaminant) #identified 266 as contaminants 
head(which(contamdf_combined$contaminant)) 

## Make presence/absence plot ====

### Make phyloseq object of presence-absence in negative controls and true samples
ps.presAbs <- transform_sample_counts(ps.All, function(abund) 1*(abund>0))
ps.presAbs.neg <- prune_samples(sample_data(ps.presAbs)$Sample_control == "Negative Control", ps.presAbs)
ps.presAbs.pos <- prune_samples(sample_data(ps.presAbs)$Sample_control == "True Sample", ps.presAbs)

# Make data.frame of prevalence in positive and negative samples
df.presAbs <- data.frame(presAbs.pos=taxa_sums(ps.presAbs.pos), presAbs.neg=taxa_sums(ps.presAbs.neg),
                         contaminant=contamdf_combined$contaminant)
plot.presAbs <- ggplot(data=df.presAbs, aes(x=presAbs.neg, y=presAbs.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence/Frequency (Negative Controls)") + ylab("Prevalence/Frequency (True Samples)")

## Save plot ====
ggplot2::ggsave(here::here("Data/03 - Decontam - Output/plot_presAbs.png"), plot.presAbs,
                height = 200, width = 200, units = "mm",
                scale = 0.5, dpi = 1000)

## Remove contaminants ====
ps.nocontam <- prune_taxa(!contamdf_combined$contaminant, ps.All) # Contains 3980 ASVs after contaminants removed

summary(taxa_sums(ps.nocontam)) 

sum(ps.nocontam@otu_table) 
sum(ps.All@otu_table) - sum(ps.nocontam@otu_table)
saveRDS(ps.nocontam, here::here("Data/03 - Decontam - Output/ps.nocontam.rds"))

## Make phyloseq object with NCs removed ====
subset_samples(ps.nocontam, Antibiotic != "Negative Control") -> ps.noneg

### Phyloseq object summary
ps.noneg <- prune_taxa(taxa_sums(ps.noneg@otu_table) > 0, ps.noneg) 

summary(taxa_sums(ps.noneg@otu_table))
sort(sample_sums(ps.noneg@otu_table)) 
