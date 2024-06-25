# Author: Sunni Patton
# Last edited: 2/1/24
# Title: Creating initial phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(phyloseq)

## Load metadata file ====
readRDS(here::here("Data/01 - Metadata - Output/sample_df_G7.rds")) -> sample_df_G7

## Make sample numbers into sample names to work with phyloseq ====
x <- sample_df_G7$Samples
rownames(sample_df_G7) <- x

## Load sequence table and taxonomy table ====
readRDS(here::here("Data/00 - Read Preprocessing - Output/seq_table_noNA_G7.rds")) -> seq_table_noNA_G7
readRDS(here::here("Data/00 - Read Preprocessing - Output/taxonomy_noNA_G7.rds")) -> taxonomy_noNA_G7

## Create phyloseq object ====
ps.All <- phyloseq(otu_table(seq_table_noNA_G7, taxa_are_rows=FALSE), sample_data(sample_df_G7), tax_table(taxonomy_noNA_G7))
saveRDS(ps.All, here::here("Data/02 - Phyloseq Preprocessing - Output/ps.All.rds"))

## Inspect phyloseq object ====
ps.All 
# Taxa distribution 
summary(taxa_sums(ps.All@otu_table)) 
# Sample read distribution
summary(sample_sums(ps.All@otu_table)) 