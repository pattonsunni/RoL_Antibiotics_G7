# Author: Sunni Patton
# Last edited: 2/1/24 
# Title: Sample processing in phyloseq
# Overview: Adding phylogenetic tree to phyloseq object and randomly subsampling to an even depth

# Note: This script involves work outside of RStudio

## Set seed ====
set.seed(123)

## Load libraries ====
library(speedyseq)
library(tidyr)
library(phytools)
library(seqateurs)
library(vegan)
library(microViz)

## Load phyloseq object (contaminants and NCs removed) ====
readRDS(here::here("Data/03 - Decontam - Output/ps.noneg.rds")) -> ps.noneg

## Add ASV column to taxonomy table ====
ps.noneg <- ps.noneg %>% mutate_tax_table(ASV = paste0("ASV", 1:3956))

## Save taxonomy table (only DNA sequence and ASV columns) ====
### Save taxonomy table as tibble
as_tibble(ps.noneg@tax_table) -> taxa_tibble
### Save only relevent columns
taxa_df <- data.frame(taxa_tibble$.otu, taxa_tibble$ASV)
### Rename columns
colnames(taxa_df) <- c("Sequence", "ASV")

readr::write_csv(taxa_df, here::here("Data/05 - Phyloseq - Output/sequenceASV.csv"))

## Change DNA sequence to ASV number ====
dna <- Biostrings::DNAStringSet(taxa_names(ps.noneg))
names(dna) <- taxa_names(ps.noneg)
ps.noneg <- merge_phyloseq(ps.noneg, dna)
taxa_names(ps.noneg) <- paste0("ASV", seq(ntaxa(ps.noneg)))

saveRDS(ps.noneg, here::here("Data/05 - Phyloseq - Output/ps.noneg.rds"))

## Load phylogenetic tree ====
# Tree created using FastTreeMP outside of Rstudio
phyloseq::read_tree(here::here("Data/treeOutput.nwk")) -> phyloTree

# Reroot the tree from midpoint 
phytools::midpoint_root(phyloTree) -> phyloTree
saveRDS(phyloTree, here::here("Data/05 - Phyloseq - Output/phyloTree.rds"))

## Add phylogenetic tree to phyloseq object ====
ps.Tree <- merge_phyloseq(ps.noneg, phyloTree)
summary(taxa_sums(ps.Tree@otu_table))
summary(sample_sums(ps.Tree@otu_table))

## Make rarefaction curve ====
as.matrix(as.data.frame(ps.Tree@otu_table)) -> data
## Supplemental Figure 3
rarecurve(data, step = 100, label = FALSE) 

## Remove samples with <1000 reads ====
sort(sample_sums(ps.Tree@otu_table))

ps.Tree <- subset_samples(ps.Tree, (Samples != "G7_T96_B3_c"))
ps.Tree <- subset_samples(ps.Tree, (Samples != "G7_T96_C2_1_b"))
ps.Tree <- subset_samples(ps.Tree, (Samples != "G7_T12_M1_2_d"))
ps.Tree <- subset_samples(ps.Tree, (Samples != "G7_T96_M2_2_b"))
ps.Tree <- subset_samples(ps.Tree, (Samples != "G7_T48_M1_1_a"))

ps.Tree <- subset_samples(ps.Tree, (Samples != "G7_T0_M2_1_a"))
ps.Tree <- subset_samples(ps.Tree, (Samples != "G7_T96_M2_3_d"))

ps.Tree <- prune_taxa(taxa_sums(ps.Tree@otu_table) > 0, ps.Tree)

summary(taxa_sums(ps.Tree@otu_table)) 

saveRDS(ps.Tree, here::here("Data/05 - Phyloseq - Output/ps.Tree.rds"))

### Prune taxa that are only present in 1 sample and with total reads < 1st quartile (56) ====
keepTaxa <- prevdt[((Prevalence >1 | TotalCounts >56)), taxaID]
ps.pruned <- prune_taxa(keepTaxa, ps.Tree)
summary(taxa_sums(ps.pruned@otu_table))

sum(ps.Tree@otu_table) - sum(ps.pruned@otu_table)

saveRDS(ps.pruned, here::here("Data/05 - Phyloseq - Output/ps.pruned.rds"))

## Randomly subsample to even depth ====
otu_table <- as.matrix(as.data.frame(ps.pruned@otu_table))
View(otu_table)
rarefied_df <- rrarefy(otu_table, sample = 5000)
sort(rowSums(rarefied_df)) 

### Rarefied ps object
rare_samData <- ps.pruned@sam_data
rare_taxTable <- ps.pruned@tax_table
rare_otuTable <- rarefied_df
rare_tree <- ps.pruned@phy_tree

ps.rare <- phyloseq(otu_table(rare_otuTable, taxa_are_rows=FALSE),sample_data(rare_samData),tax_table(rare_taxTable), phy_tree(rare_tree))

ps.rare <- prune_taxa(taxa_sums(ps.rare@otu_table) > 0, ps.rare)

summary(taxa_sums(ps.rare))

saveRDS(ps.rare, here::here("Data/05 - Phyloseq - Output/ps.rare.rds"))