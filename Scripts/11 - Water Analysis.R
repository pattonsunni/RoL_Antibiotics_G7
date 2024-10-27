# Author: Denise Silva (Patton et al. 2024)
# Title: Water Analysis
# Date: 10/27/24 (last edited)

# Overview: 
# Loaded in all fasta files, ran all quality control steps, and inferred sample sequences and taxonomy. Created sequence table and taxonomy table that will be used later when creating a phyloseq object.
# Built a metadata file based on original file name structure to use in phyloseq object
# Created initial phyloseq object
# Removed contaminants

## Load libraries ====
library(here)
library(dada2)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(decontam)
library(ggplot2)
library(seqateurs)
library(speedyseq)
library(taxreturn)


## Set seed for reproducibility ====
set.seed(123)

## Read Preprocessing ====

# Define file paths
path_G7 <- here::here("./processing_g7/raw_sequences")
fnFs_G7 <- sort(list.files(path_G7, pattern = "R1_001.fastq.gz", 
                           full.names = TRUE))
fnRs_G7 <- sort(list.files(path_G7, pattern = "R2_001.fastq.gz", 
                           full.names = TRUE))

# Extract sample names from files
# Assume name follows: lane1-s0001-index--Findenx-Rindex-SAMPLENAME_XXX.fastq
sampleNames_G7 <- 
  sapply(strsplit(basename(fnFs_G7), "-"), `[`,7) 

# Removes all beginning information, but still left with information at the end
# (S#_R1_001.fastq.gz)
sampleNames_G7 <- 
  sapply(strsplit(basename(sampleNames_G7), "\\."), `[`,1) 

# Removes fastq.gz at the end, but still left with the _S*_R1_001
sampleNames_G7 <- gsub("_S\\d+\\d+\\d+\\d+_R1_001", "", sampleNames_G7) 

# Remove "_stvx" from the samples
sampleNames_G7 <- str_replace(sampleNames_G7, "_stvx", "")

# Set paths for quality-filtered files
filtFs_G7 <- file.path(path_G7, "filterAndTrim", paste0(sampleNames_G7,
                                                        "_F_filt.fastq.gz"))
filtRs_G7 <- file.path(path_G7, "filterAndTrim", paste0(sampleNames_G7,
                                                        "_R_filt.fastq.gz"))
# Perform Quality filtering and trimming
out_G7 <- filterAndTrim(fnFs_G7, filtFs_G7, fnRs_G7, filtRs_G7, truncLen = c(245, 230),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=FALSE)

saveRDS(out_G7, here::here("./processing_g7/outputfiles_g7/00_read_processing/out_G7.rds"))

# Assess number of reads loss
sum(out_G7[,1])-sum(out_G7[,2]) 

# Filter ratio
filter_fun = data.frame(out_G7)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio 

## Learn errors
errF_G7 <- learnErrors(filtFs_G7,multithread = FALSE, nbases = 5e8) 
plotErrors(errF_G7, nominalQ=TRUE)

saveRDS(errF_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/errF_G7.rds"))

errR_G7<-learnErrors(filtRs_G7,multithread = FALSE, nbases = 5e8) 
plotErrors(errR_G7, nominalQ=TRUE)

saveRDS(errR_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/errR_G7.rds"))

# Ensure sample naming is consistent
sampleNames_G7 <- names(filtFs_G7)
sampleNames_G7 <-names(filtRs_G7)

# Infer sample sequence (dada2)
dadaForward_G7 <- 
  dada(filtFs_G7, 
       err=errF_G7, 
       multithread=FALSE)
saveRDS(dadaForward_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/dadaForward_G7.rds"))

dadaReverse_G7 <- 
  dada(filtRs_G7, err=errR_G7, multithread=FALSE)
saveRDS(dadaReverse_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/dadaReverse_G7.rds"))

# Merge reads into contigs
contigs_G7 <-
  mergePairs(dadaForward_G7, filtFs_G7, dadaReverse_G7, filtRs_G7)

saveRDS(contigs_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/contigs_G7.rds"))

# Sequence table and filter by contig length
seq_table_G7 <-
  makeSequenceTable(contigs_G7) 
dim(seq_table_G7) 

table(nchar(getSequences(seq_table_G7))) 

seq_table_G7<-
  seq_table_G7[,nchar(colnames(seq_table_G7)) %in% 253:254]
table(nchar(getSequences(seq_table_G7)))
dim(seq_table_G7) 
sum(seq_table_G7) 

saveRDS(seq_table_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/seq_table_G7.rds"))

# Remove chimeras
seq_table_nochimeri_G7 <- 
  removeBimeraDenovo(seq_table_G7, method="consensus", multithread=TRUE, verbose=TRUE) 

dim(seq_table_nochimeri_G7) 

sum(seq_table_G7) - sum(seq_table_nochimeri_G7) 

saveRDS(seq_table_nochimeri_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/seq_table_nochimeri_G7.rds"))

# Assign taxonomy
taxa_G7 <- assignTaxonomy(seq_table_nochimeri_G7, here::here("./Training_Naive_bayes/silva_nr99_v138.1_train_set.fa.gz"), multithread=FALSE)
taxa_G7 <- addSpecies(taxa_G7, here::here("./Training_Naive_bayes/silva_species_assignment_v138.1.fa.gz"))

saveRDS(taxa_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/taxa_G7.rds")) 

# Remove chloroplast and mitochondrial sequences
is.chloroplast_G7 <- taxa_G7[,"Order"] %in% "Chloroplast"
seq_table_nochloro_G7 <- seq_table_nochimeri_G7[,!is.chloroplast_G7]
dim(seq_table_nochloro_G7) 
sum(seq_table_nochimeri_G7) - sum(seq_table_nochloro_G7) 
saveRDS(seq_table_nochloro_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/seq_table_nochloro_G7.rds"))

# New taxonomy table (no chloroplast)
taxonomy_nochloro_G7 <- taxa_G7[!is.chloroplast_G7,]
dim(taxonomy_nochloro_G7) 
View(taxonomy_nochloro_G7)
saveRDS(taxonomy_nochloro_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/taxonomy_nochlor_G7.rds"))

is.mitochondria_G7 <- 
  taxonomy_nochloro_G7[,"Family"] %in% "Mitochondria"
seq_table_nomito_G7 <- 
  seq_table_nochloro_G7[,!is.mitochondria_G7]
dim(seq_table_nomito_G7) 
sum(seq_table_nochloro_G7) - sum(seq_table_nomito_G7) 
saveRDS(seq_table_nomito_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/seq_table_nomito_G7.rds"))

# New taxonomy table (no mitochondria and chloroplast)
taxonomy_nomito_G7 <-
  taxonomy_nochloro_G7[!is.mitochondria_G7,]
dim(taxonomy_nomito_G7) 
View(taxonomy_nomito_G7)
saveRDS(taxonomy_nomito_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/taxonomy_nomito_G7.rds"))

# Remove sequences not annotated beyond the Kingdom level
is.NA_G7 <- 
  taxonomy_nomito_G7[,"Phylum"] %in% NA
seq_table_noNA_G7 <- seq_table_nomito_G7[,!is.NA_G7]
dim(seq_table_noNA_G7) 

sum(seq_table_nomito_G7) - sum(seq_table_noNA_G7) 
saveRDS(seq_table_noNA_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/seq_table_noNA_G7.rds"))

taxonomy_noNA_G7<- taxonomy_nomito_G7[!is.NA_G7,]
dim(taxonomy_noNA_G7) 
saveRDS(taxonomy_noNA_G7, 
        here::here("./processing_g7/outputfiles_g7/00_read_processing/taxonomy_noNA_G7.rds"))

# Assess total number of reads removed and ASV frequency
sum(seq_table_noNA_G7) 
dim(seq_table_noNA_G7)  

summary(colSums(seq_table_noNA_G7)) 

summary(rowSums(seq_table_noNA_G7)) 

# See which sample had that low read number
sort(rowSums(seq_table_noNA_G7))

# Tracking reads removed at each step
getN <- function(x) sum(getUniques(x))
track_G7 <- cbind(out_G7, 
                  sapply(dadaForward_G7, getN), 
                  sapply(dadaReverse_G7, getN), 
                  sapply(contigs_G7, getN), 
                  rowSums(seq_table_nochimeri_G7),
                  rowSums(seq_table_nochloro_G7), 
                  rowSums(seq_table_nomito_G7), 
                  rowSums(seq_table_noNA_G7))
colnames(track_G7) <- c("input", 
                        "filtered", 
                        "denoisedF", 
                        "denoisedR", "merged",
                        "nonchim", 
                        "nochloro", 
                        "nomito", 
                        "noNA.phylum")

# Update sample names based on track_G7 rows 
sampleNames_G7_new <- sapply(strsplit(basename(rownames(track_G7)), "-"), `[`,7) 
sampleNames_G7_new <- sapply(strsplit(basename(sampleNames_G7_new), "\\."), `[`,1) 
sampleNames_G7_new <- gsub("_S\\d+\\d?\\d?\\d?_R1_001", "", sampleNames_G7_new) 
sampleNames_G7_new

# Save sampleNames_G7_new 
saveRDS(sampleNames_G7_new, here::here("./processing_g7/outputfiles_g7/00_read_processing/sampleNames.rds"))

# Assign new row names (shorter sample name) using sampleNames_G7_new
rownames(track_G7) <- sampleNames_G7_new

write.csv(track_G7, here::here("./processing_g7/outputfiles_g7/00_read_processing/track_G7.csv")) 

readr::read_csv(here::here("./processing_g7/outputfiles_g7/00_read_processing/track_G7.csv"), col_names = TRUE) -> track_G7 # Considers the sample name column a true column, but has weird name (...1)

# Rename column
x <- c(Samples = "...1")
x <- c(x, new = "Samples")

track_G7 <- 
  rename(track_G7, any_of(x))

# Save new CSV with correct Samples column
readr::write_csv(track_G7, 
                 here::here("./processing_g7/outputfiles_g7/00_read_processing/track_G7.csv"))

# Library size data frame
librarySize_df <- 
  data.frame(track_G7$Samples, 
             track_G7$input, 
             track_G7$noNA.phylum)

# Rename columns
colnames(librarySize_df) <- c("Samples", "PreQC", "PostQC")

readr::write_csv(librarySize_df, 
                 here::here("./processing_g7/outputfiles_g7/00_read_processing/librarySize.csv"))



## Metadata ====
sampleNames_G7_new <- 
  readRDS("./processing_g7/outputfiles_g7/00_read_processing/sampleNames.rds")

# Remove "_stvx" from the samples
sampleNames_G7_new <- 
  str_replace(sampleNames_G7_new, "_stvx", "")

# Extract information from sampleNames_G7_new
# Removed treatment, tank, and replicate
time_G7 <- 
  as.character(sapply(strsplit(sampleNames_G7_new, "G7_T"), `[`,2)) 

time_G7 <- 
  as.character(sapply(strsplit(time_G7, "_"), `[`, 1)) 

treatment_G7 <- 
  as.character(sapply(strsplit(sampleNames_G7_new, "G7_T\\d+\\_"), `[`, 2)) 

treatment_G7 <-
  as.character(sapply(strsplit(treatment_G7, "_"), `[`, 1)) 

treatment.rep_G7 <-
  as.character(sapply(strsplit(sampleNames_G7_new, "G7_T\\d+\\_"), `[`, 2))

treatment.rep_G7 <- 
  as.character(sapply(strsplit(treatment.rep_G7, "_[a-d]"), `[`, 1))

# Make dataframe
sample_water_df_G7 <- data.frame(Samples = sampleNames_G7_new,
                                 Time = time_G7, 
                                 Treatment = treatment_G7, 
                                 TreatmentTank = treatment.rep_G7)

# Change Treatment's name
sample_water_df_G7 <- sample_water_df_G7 %>%
  mutate(Treatment = if_else(Treatment %in% c("B1", 
                                              "B2", 
                                              "B3", 
                                              "B4", 
                                              "B5", 
                                              "B6"), 
                             "Blank", Treatment))

sample_water_df_G7 <- sample_water_df_G7 %>%
  mutate(Treatment = recode(Treatment,
                            "M1" = "Mixture Low",
                            "M2" = "Mixture High",
                            "A1" = "Ampicillin Low",
                            "A2" = "Ampicillin High",
                            "S1" = "Streptomycin Low",
                            "S2" = "Streptomycin High",
                            "C1" = "Ciprofloxacin Low",
                            "C2" = "Ciprofloxacin High"))

# Add Antibiotic column
# Creating a new column and using negative control as a placeholder.
sample_water_df_G7$Antibiotic <- "NA" 

# Add correct names-values
sample_water_df_G7 <- sample_water_df_G7 %>%
  mutate(Antibiotic = case_when(
    Treatment == "Blank" ~ "Blank",
    Treatment %in% c("Mixture Low", "Mixture High") ~ "Mixture",
    Treatment %in% c("Ampicillin Low", "Ampicillin High") ~ "Ampicillin",
    Treatment %in% c("Streptomycin Low", "Streptomycin High") ~ "Streptomycin",
    Treatment %in% c("Ciprofloxacin Low", "Ciprofloxacin High") ~ "Ciprofloxacin"))

# Add Dose column
# Using Blank as placeholder
sample_water_df_G7$Dose <- "Blank" 

# Add information to Dose column
sample_water_df_G7 <- sample_water_df_G7 %>%
  mutate(Dose = case_when(
    Treatment %in% c("Mixture Low", 
                     "Ampicillin Low", 
                     "Streptomycin Low", 
                     "Ciprofloxacin Low") ~ "Low",  
    TRUE ~ Dose
  ))

sample_water_df_G7 <- sample_water_df_G7 %>%
  mutate(Dose = case_when(
    Treatment %in% c("Mixture High",
                     "Ampicillin High",
                     "Streptomycin High",
                     "Ciprofloxacin High") ~ "High",
    TRUE ~ Dose
  ))

# Add Pretreatment column(s)
# All T0 groups are actually pretreatment 
sample_water_df_G7$Pretreatment <- sample_water_df_G7$Treatment

# Add information to Pretreatment column
sample_water_df_G7$Pretreatment[sample_water_df_G7$Time == "0"] <- "Blank"

# Add PretreatmentAbx column
sample_water_df_G7$PretreatmentAbx <- sample_water_df_G7$Antibiotic

# Add information to PretreatmentAbx column
sample_water_df_G7$PretreatmentAbx[sample_water_df_G7$Time == "0"] <- "Blank"

# Add PretreamentDose column
sample_water_df_G7$PretreatmentDose <- sample_water_df_G7$Dose

# Add information to PretreatmentAbx column
sample_water_df_G7$PretreatmentDose[sample_water_df_G7$Time == "0"] <- "Blank"

# Add Tank column
# Designating separate tanks
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "B1"] <-"1"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "B2"] <-"2"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "B3"] <-"3"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "B4"] <-"4"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "B5"] <-"5"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "B6"] <-"6"

sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "A1_1"] <-"7"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "A1_2"] <-"8"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "A1_3"] <-"9"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "A2_1"] <-"10"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "A2_2"] <-"11"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "A2_3"] <-"12"

sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "S1_1"] <-"13"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "S1_2"] <-"14"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "S1_3"] <-"15"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "S2_1"] <-"16"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "S2_2"] <-"17"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "S2_3"] <-"18"

sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "C1_1"] <-"19"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "C1_2"] <-"20"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "C1_3"] <-"21"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "C2_1"] <-"22"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "C2_2"] <-"23"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "C2_3"] <-"24"

sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "M1_1"] <-"25"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "M1_2"] <-"26"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "M1_3"] <-"27"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "M2_1"] <-"28"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "M2_2"] <-"29"
sample_water_df_G7$Tank[sample_water_df_G7$TreatmentTank == "M2_3"] <-"30"

sample_water_df_G7 <- sample_water_df_G7 %>%
  mutate(across(everything(), ~ if_else(is.na(.), 
                                        "Negative Control", as.character(.))))

# Additional information for decontam and library size plotting)
# Add Sample_control column
sample_water_df_G7$Sample_control <- "True Sample"

# Add information to Sample_control column
sample_water_df_G7$Sample_control[sample_water_df_G7$Antibiotic == "Negative Control"] <- "Negative Control"

# Add DNA quantification column 
quant_G7 <- 
  readr::read_csv(here::here("./processing_g7/quant_G7.csv")) 

# Join the two dataframes 
sample_water_df_G7 <- sample_water_df_G7 %>%
  left_join(quant_G7, by = "Samples")

# Add LibrarySize columns (pre and post QC)
librarySize_water_df <- 
  readr::read_csv(here::here("./processing_g7/outputfiles_g7/00_read_processing/librarySize.csv"))  

# Join the dataframes 
sample_water_df_G7 <- sample_water_df_G7 %>%     
  left_join(librarySize_water_df, by = "Samples")

# Save metadata file
saveRDS(sample_water_df_G7, here::here("./processing_g7/outputfiles_g7/01_metadata/sample_water_df_G7.rds"))
readr::write_csv(sample_water_df_G7, here::here("./processing_g7/outputfiles_g7/01_metadata/sample_water_df_G7.csv"))


## Phyloseq Preprocessing ====
# Load metadata file
sample_water_df_G7 <- 
  readRDS("./processing_g7/outputfiles_g7/01_metadata/sample_df_G7.rds")

# Make sample numbers into sample names
rownames(sample_water_df_G7) <- sample_water_df_G7$Samples

# Load sequence table and taxonomy table
seq_table_water_noNA_G7 <- 
  readRDS("./processing_g7/outputfiles_g7/00_read_processing/seq_table_noNA_G7.rds")

rownames(seq_table_water_noNA_G7) <- 
  rownames(seq_table_water_noNA_G7) %>% 
  str_replace("_stvx_F_filt.fastq.gz|_F_filt.fastq.gz", "")

taxonomy_water_noNA_G7 <- 
  readRDS("./processing_g7/outputfiles_g7/00_read_processing/taxonomy_noNA_G7.rds")

# Create phyloseq object
ps.All <- phyloseq(otu_table(seq_table_water_noNA_G7,
                             taxa_are_rows=FALSE),
                   sample_data(sample_water_df_G7), 
                   tax_table(taxonomy_water_noNA_G7))

saveRDS(ps.All, ("./processing_g7/outputfiles_g7/02_phylo_processing/ps.All.rds"))

# Taxa distribution 
summary(taxa_sums(ps.All@otu_table)) 

# Sample read distribution
summary(sample_sums(ps.All@otu_table)) 

## Decontam ====
# Plot library size by sample
df <- as.data.frame(sample_data(ps.All)) 
df$LibrarySize <- sample_sums(ps.All)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot_decontam <-
  ggplot(
    data=df, 
    aes(x=Index,
        y=LibrarySize, 
        color=Sample_control)) + geom_point()

# Save plot
ggplot2::ggsave(("./processing_g7/outputfiles_g7/data_003_decontam/plot_decontam.png"),
                plot_decontam,
                height = 200, width = 350, units = "mm",
                scale = 0.5, dpi = 1000)

# Run decontam (combined)
# Add a column called is.neg (TRUE/ FALSE)
sample_data(ps.All)$is.neg <- 
  sample_data(ps.All)$Sample_control == "Negative Control" 

contamdf_combined <- isContaminant(ps.All, method="combined", 
                                   neg="is.neg", 
                                   conc="quantification", 
                                   threshold=0.5)

head(which(contamdf_combined$contaminant)) 

# Make presence/absence plot
# From https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Make phyloseq object of presence-absence in negative controls and true samples
ps.presAbs <- transform_sample_counts(ps.All, function(abund) 1*(abund>0))
ps.presAbs.neg <- prune_samples(sample_data(ps.presAbs)$Sample_control == "Negative Control",
                                ps.presAbs)
ps.presAbs.pos <- prune_samples(sample_data(ps.presAbs)$Sample_control == "True Sample", 
                                ps.presAbs)

# Make data.frame of prevalence in positive and negative samples
df.presAbs <- data.frame(presAbs.pos=taxa_sums(ps.presAbs.pos), presAbs.neg=taxa_sums(ps.presAbs.neg),
                         contaminant=contamdf_combined$contaminant)
plot.presAbs <- ggplot(data=df.presAbs,
                       aes(x=presAbs.neg, 
                           y=presAbs.pos, 
                           color=contaminant)) + 
  geom_point() +
  xlab("Prevalence/Frequency (Negative Controls)") +
  ylab("Prevalence/Frequency (True Samples)")

## Save plot
ggplot2::ggsave(here::here("./processing_g7/outputfiles_g7/data_003_decontam/plot_presAbs.png"), 
                plot.presAbs,
                height = 200, 
                width = 200, 
                units = "mm",
                scale = 0.5,
                dpi = 1000)

## Remove contaminants
ps.nocontam <- prune_taxa(!contamdf_combined$contaminant, ps.All) 

summary(taxa_sums(ps.nocontam)) 
sum(ps.nocontam@otu_table) 
sum(ps.All@otu_table) - sum(ps.nocontam@otu_table) 
saveRDS(ps.nocontam, here::here("./processing_g7/outputfiles_g7/data_003_decontam/ps.nocontam.rds"))

ps.nocontam <- 
  readr::read_rds("./processing_g7/outputfiles_g7/data_003_decontam/ps.nocontam.rds")

## Make phyloseq object with Negative control samples removed
ps.noneg <- 
  subset_samples(ps.nocontam, Antibiotic != "Negative Control") 

# Phyloseq object summary
ps.noneg <- 
  prune_taxa(taxa_sums(ps.noneg@otu_table) > 0, ps.noneg) 

summary(taxa_sums(ps.noneg@otu_table))

sort(sample_sums(ps.noneg@otu_table)) 

saveRDS(ps.noneg, here::here("./processing_g7/outputfiles_g7/data_003_decontam/ps.noneg.rds"))