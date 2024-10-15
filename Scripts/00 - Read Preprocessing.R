# Author: Sunni Patton
# Last edited: 01/30/24
# Title: Read preprocessing
# Overview: Loading raw data and preprocessing data following the DADA2 pipeline. Taxonomy assigned using Silva 138.1

## Load libraries ====
library(here)
library(dada2)
library(dplyr)

## Set seed for reproducibility ====
set.seed(123)

## Pre-Processing ====

## Set path to raw files
### These are raw reads that have been run through cutadapt script (on terminal) to remove primers
path_G7 <- here::here("Data/Raw Reads")

### Set path for forward and reverse reads
fnFs_G7 <- sort(list.files(path_G7, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs_G7 <- sort(list.files(path_G7, pattern = "R2_001.fastq.gz", full.names = TRUE))

## Extract sample names from files
### Assume name follows: lane1-s0001-index--Findenx-Rindex-SAMPLENAME_XXX.fastq
sampleNames_G7 <- sapply(strsplit(basename(fnFs_G7), "-"), `[`,7) 
sampleNames_G7 <- sapply(strsplit(basename(sampleNames_G7), "\\."), `[`,1) 
sampleNames_G7 <- gsub("_S\\d+\\d?\\d?\\d?_R1_001", "", sampleNames_G7) 

## Set file destination for files after quality filtering
filtFs_G7 <- file.path(path_G7, "filterAndTrim", paste0(sampleNames_G7, "_F_filt.fastq.gz"))
filtRs_G7 <- file.path(path_G7, "filterAndTrim", paste0(sampleNames_G7, "_R_filt.fastq.gz"))

## Filter and Trim ====
out_G7 <- filterAndTrim(fnFs_G7, filtFs_G7, fnRs_G7, filtRs_G7, truncLen = c(245, 230),
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=FALSE)
## Save output
saveRDS(out_G7, here::here("/Data/00 - Read Preprocessing - Output/out_G7.rds"))

## Assess number of reads lost
sum(out_G7[,1])-sum(out_G7[,2]) 

filter_fun = data.frame(out_G7)
filter_fun
ratio = sum(filter_fun$reads.out)/sum(filter_fun$reads.in)
ratio

## Learn errors and infer sample sequence ====
errF_G7 <- learnErrors(filtFs_G7,multithread = FALSE, nbases = 5e8)  
#504,570,395 total bases in 2059471 reads from 83 samples will be used for learning the error rates.
saveRDS(errF_G7, here::here("Data/00 - Read Preprocessing - Output/errF_G7.rds"))

errR_G7<-learnErrors(filtRs_G7,multithread = FALSE, nbases = 5e8) 
#500,388,000  total bases in 2,175,600 reads from 88 samples samples will be used for learning the error rates.
saveRDS(errR_G7, here::here("Data/00 - Read Preprocessing - Output/errR_G7.rds"))

## Ensure sample naming is consistent
names(filtFs_G7)<-sampleNames_G7
names(filtRs_G7)<-sampleNames_G7

## Infer sample sequence
dadaForward_G7 <- dada(filtFs_G7, err=errF_G7, multithread=FALSE)
saveRDS(dadaForward_G7, here::here("/Data/00 - Read Preprocessing - Output/dadaForward_G7.rds"))

dadaReverse_G7 <- dada(filtRs_G7, err=errR_G7, multithread=FALSE)
saveRDS(dadaReverse_G7, here::here("/Data/00 - Read Preprocessing - Output/dadaReverse_G7.rds"))

## Create contigs and sequence table ====
contigs_G7 <- mergePairs(dadaForward_G7, filtFs_G7, dadaReverse_G7, filtRs_G7)
saveRDS(contigs_G7, here::here("/Data/00 - Read Preprocessing - Output/contigs_G7.rds"))

## Make sequence table and visualize contig length and frequency
seq_table_G7 <- makeSequenceTable(contigs_G7) 
dim(seq_table_G7)

table(nchar(getSequences(seq_table_G7))) 

### Keep contigs within desired size range
seq_table_G7<-seq_table_G7[,nchar(colnames(seq_table_G7)) %in% 251:255]
table(nchar(getSequences(seq_table_G7)))
dim(seq_table_G7) 
sum(seq_table_G7) 

## Save output
saveRDS(seq_table_G7, here::here("/Data/00 - Read Preprocessing - Output/seq_table_G7.rds"))

## Remove chimeras ====
seq_table_nochimeri_G7 <- removeBimeraDenovo(seq_table_G7, method="consensus", multithread=TRUE, verbose=TRUE) 
dim(seq_table_nochimeri_G7) 
sum(seq_table_G7) - sum(seq_table_nochimeri_G7)

saveRDS(seq_table_nochimeri_G7, here::here("/Data/00 - Read Preprocessing - Output/seq_table_nochimeri_G7.rds"))

## Assign taxonomy ====
taxa_G7 <- assignTaxonomy(seq_table_nochimeri_G7, here::here("Data/silva_nr99_v138.1_train_set.fa.gz"), multithread=FALSE)
taxa_G7 <- addSpecies(taxa_G7, here::here("Data/silva_species_assignment_v138.1.fa.gz"))

saveRDS(taxa_G7, here::here("Data/00 - Read Preprocessing - Output/taxa_G7.rds")) 

## Remove off-target sequences ====
### New sequence table (no chloroplast)
is.chloroplast_G7 <- taxa_G7[,"Order"] %in% "Chloroplast"
seq_table_nochloro_G7 <- seq_table_nochimeri_G7[,!is.chloroplast_G7]
dim(seq_table_nochloro_G7) 
sum(seq_table_nochimeri_G7) - sum(seq_table_nochloro_G7) 
saveRDS(seq_table_nochloro_G7, here::here("/Data/00 - Read Preprocessing - Output/seq_table_nochloro_G7.rds"))

### New taxonomy table (no chloroplast)
taxonomy_nochloro_G7 <- taxa_G7[!is.chloroplast_G7,]
dim(taxonomy_nochloro_G7)
View(taxonomy_nochloro_G7)
saveRDS(taxonomy_nochloro_G7, here::here("Data/00 - Read Preprocessing - Output/taxonomy_nochlor_G7.rds"))

### New sequence table (no mitochondria or chloroplast)
is.mitochondria_G7 <- taxonomy_nochloro_G7[,"Family"] %in% "Mitochondria"
seq_table_nomito_G7 <- seq_table_nochloro_G7[,!is.mitochondria_G7]
dim(seq_table_nomito_G7)
sum(seq_table_nochloro_G7) - sum(seq_table_nomito_G7) 
saveRDS(seq_table_nomito_G7, here::here("Data/00 - Read Preprocessing - Output/seq_table_nomito_G7.rds"))

### New taxonomy table (no mitochondria or chloroplast)
taxonomy_nomito_G7 <- taxonomy_nochloro_G7[!is.mitochondria_G7,]
dim(taxonomy_nomito_G7)
View(taxonomy_nomito_G7)
saveRDS(taxonomy_nomito_G7, here::here("Data/00 - Read Preprocessing - Output/taxonomy_nomito_G7.rds"))

## Remove sequences not annotated beyond the Kingdom level
is.NA_G7 <- taxonomy_nomito_G7[,"Phylum"] %in% NA
seq_table_noNA_G7 <- seq_table_nomito_G7[,!is.NA_G7]
dim(seq_table_noNA_G7) 

sum(seq_table_nomito_G7) - sum(seq_table_noNA_G7) 
saveRDS(seq_table_noNA_G7, here::here("/Data/00 - Read Preprocessing - Output/seq_table_noNA_G7.rds"))

taxonomy_noNA_G7<- taxonomy_nomito_G7[!is.NA_G7,]
dim(taxonomy_noNA_G7)
saveRDS(taxonomy_noNA_G7, here::here("/Data/00 - Read Preprocessing - Output/taxonomy_noNA_G7.rds"))

## Assess total number of reads removed and ASV frequency ====
sum(seq_table_noNA_G7) # 6,425,441 reads left after quality control
dim(seq_table_noNA_G7) # 4246 contigs left after quality control

summary(colSums(seq_table_noNA_G7)) 
summary(rowSums(seq_table_noNA_G7)) # Minimum reads in a sample is 136 

sort(rowSums(seq_table_noNA_G7)) # Sample with 136 reads is a negative control 

## Tracking reads removed at each step ====
getN <- function(x) sum(getUniques(x))
track_G7 <- cbind(out_G7, sapply(dadaForward_G7, getN), sapply(dadaReverse_G7, getN), sapply(contigs_G7, getN), rowSums(seq_table_nochimeri_G7), rowSums(seq_table_nochloro_G7), rowSums(seq_table_nomito_G7), rowSums(seq_table_noNA_G7))
colnames(track_G7) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "nochloro", "nomito", "noNA.phylum")

### Update sample names based on track_G7 rows (simply forcing column names to be sampleNames_G7 incorrectly assigns the names)
sampleNames_G7_new <- sapply(strsplit(basename(rownames(track_G7)), "-"), `[`,7) # Removes all beginning information, but still left with information at the end (S#_R1_001.fastq.gz)
sampleNames_G7_new <- sapply(strsplit(basename(sampleNames_G7_new), "\\."), `[`,1) # Removes fastq.gz at the end, but still left with the _S*_R1_001
sampleNames_G7_new <- gsub("_S\\d+\\d?\\d?\\d?_R1_001", "", sampleNames_G7_new) 
sampleNames_G7_new

### Save sampleNames_G7_new (will need later for metadata)
saveRDS(sampleNames_G7_new, here::here("Data/00 - Read Preprocessing - Output/sampleNames.rds"))

### Assign new row names (shorter sample name) using sampleNames_G7_new
rownames(track_G7) <- sampleNames_G7_new

write.csv(track_G7, here::here("Data/00 - Read Preprocessing - Output/track_G7.csv")) # Saves column with sample names but that column doesn't have a column name

### Read CSV back in
readr::read_csv(here::here("Data/00 - Read Preprocessing - Output/track_G7.csv"), col_names = TRUE) -> track_G7 # Considers the sample name column a true column, but has weird name (...1)

### Rename column
x <- c(Samples = "...1")
x <- c(x, new = "Samples")

track_G7 <- rename(track_G7, any_of(x))

### Save new CSV with correct Samples column
readr::write_csv(track_G7, here::here("Data/00 - Read Preprocessing - Output/track_G7.csv"))

### Save another file with sample names, library size (preQC), and library size (postQC)
librarySize_df <- data.frame(track_G7$Samples, track_G7$input, track_G7$noNA.phylum)
### Rename columns
colnames(librarySize_df) <- c("Samples", "PreQC", "PostQC")
readr::write_csv(librarySize_df, here::here("Data/00 - Read Preprocessing - Output/librarySize.csv"))
