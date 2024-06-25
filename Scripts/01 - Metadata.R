# Author: Sunni Patton
# Last edited: 2/1/24
# Title: Creating metadata file
# Overview: Built a metadata file based on original file name structure to use in phyloseq object

## Set seed ====
set.seed(123)

## Load libraries ====
library(ggplot2)
library(dplyr)
library(rstatix)
library(broom)

## Extract sample name data ====
### Read in sampleNames_G7_new if not already loaded
readRDS(here::here("Data/00 - Read Preprocessing - Output/sampleNames.rds")) -> sampleNames_G7_new

### Extract information from sampleNames_G7_new
time_G7 <- as.character(sapply(strsplit(sampleNames_G7_new, "_T"), `[`,2))
time_G7 <- as.character(sapply(strsplit(time_G7, "_"), `[`, 1)) 

treatment_G7 <- as.character(sapply(strsplit(sampleNames_G7_new, "G7_T\\d+\\_"), `[`, 2)) 
treatment_G7 <- as.character(sapply(strsplit(treatment_G7, "_"), `[`, 1)) 

fragment_G7 <- as.character(sapply(strsplit(sampleNames_G7_new, "G7_T\\d+\\_"), `[`, 2))
fragment_G7 <- as.character(sapply(strsplit(fragment_G7, "[A-Z]\\d+\\_\\d?\\_?"), `[`, 2))

treatment.rep_G7 <- as.character(sapply(strsplit(sampleNames_G7_new, "G7_T\\d+\\_"), `[`, 2))
treatment.rep_G7 <- as.character(sapply(strsplit(treatment.rep_G7, "_[a-d]"), `[`, 1))


## Create metadata dataframe ====
sample_df_G7 <- data.frame(Samples = sampleNames_G7_new, Time = time_G7, Treatment = treatment_G7, Fragment = fragment_G7, TreatmentTank = treatment.rep_G7)

## Make Treatment long 
sample_df_G7$Treatment[sample_df_G7$Treatment == "B1" | 
                         sample_df_G7$Treatment == "B2" | 
                         sample_df_G7$Treatment == "B3" | 
                         sample_df_G7$Treatment == "B4" | 
                         sample_df_G7$Treatment == "B5" | 
                         sample_df_G7$Treatment == "B6"] <- "Blank"

sample_df_G7$Treatment[sample_df_G7$Treatment == "M1"] <-"Mixture Low"
sample_df_G7$Treatment[sample_df_G7$Treatment == "M2"] <-"Mixture High"
sample_df_G7$Treatment[sample_df_G7$Treatment == "A1"] <-"Ampicillin Low"
sample_df_G7$Treatment[sample_df_G7$Treatment == "A2"] <-"Ampicillin High"
sample_df_G7$Treatment[sample_df_G7$Treatment == "S1"] <-"Streptomycin Low"
sample_df_G7$Treatment[sample_df_G7$Treatment == "S2"] <-"Streptomycin High"
sample_df_G7$Treatment[sample_df_G7$Treatment == "C1"] <-"Ciprofloxacin Low"
sample_df_G7$Treatment[sample_df_G7$Treatment == "C2"] <-"Ciprofloxacin High"

## Add Antibiotic column 
sample_df_G7$Antibiotic <- "Negative Control" # Creating a new column in our dataframe and using negative control as a placeholder

### Add correct values
sample_df_G7$Antibiotic[sample_df_G7$Treatment == "Blank"] <- "Blank"

sample_df_G7$Antibiotic[sample_df_G7$Treatment == "Mixture Low" |
                          sample_df_G7$Treatment == "Mixture High"] <- "Mixture"

sample_df_G7$Antibiotic[sample_df_G7$Treatment == "Ampicillin Low" |
                          sample_df_G7$Treatment == "Ampicillin High"] <- "Ampicillin"

sample_df_G7$Antibiotic[sample_df_G7$Treatment == "Streptomycin Low" |
                          sample_df_G7$Treatment == "Streptomycin High"] <- "Streptomycin"

sample_df_G7$Antibiotic[sample_df_G7$Treatment == "Ciprofloxacin Low" |
                          sample_df_G7$Treatment == "Ciprofloxacin High"] <- "Ciprofloxacin"


## Add Dose column
sample_df_G7$Dose <- "Blank" # Using Blank as placeholder

### Add information to Dose column
sample_df_G7$Dose[sample_df_G7$Treatment == "Mixture Low" |
                    sample_df_G7$Treatment == "Ampicillin Low" |
                    sample_df_G7$Treatment == "Streptomycin Low" |
                    sample_df_G7$Treatment == "Ciprofloxacin Low"] <- "Low"

sample_df_G7$Dose[sample_df_G7$Treatment == "Mixture High" |
                    sample_df_G7$Treatment == "Ampicillin High" |
                    sample_df_G7$Treatment == "Streptomycin High" |
                    sample_df_G7$Treatment == "Ciprofloxacin High"] <- "High"


## Add Pretreatment column(s)
### Adding these to indicate that all T0 groups are actually pretreatment 

### Add Pretreatment column
sample_df_G7$Pretreatment <- sample_df_G7$Treatment
### Add information to Pretreatment column
sample_df_G7$Pretreatment[sample_df_G7$Time == "0"] <- "Blank"

### Add PretreatmentAbx column
sample_df_G7$PretreatmentAbx <- sample_df_G7$Antibiotic
### Add information to PretreatmentAbx column
sample_df_G7$PretreatmentAbx[sample_df_G7$Time == "0"] <- "Blank"

### Add PretreamentDose column
sample_df_G7$PretreatmentDose <- sample_df_G7$Dose
### Add information to PretreatmentAbx column
sample_df_G7$PretreatmentDose[sample_df_G7$Time == "0"] <- "Blank"

## Add Tank column 
### These numbers are arbitrarily assigned to designate separate tanks

sample_df_G7$Tank[sample_df_G7$TreatmentTank == "B1"] <-"1"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "B2"] <-"2"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "B3"] <-"3"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "B4"] <-"4"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "B5"] <-"5"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "B6"] <-"6"

sample_df_G7$Tank[sample_df_G7$TreatmentTank == "A1_1"] <-"7"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "A1_2"] <-"8"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "A1_3"] <-"9"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "A2_1"] <-"10"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "A2_2"] <-"11"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "A2_3"] <-"12"

sample_df_G7$Tank[sample_df_G7$TreatmentTank == "S1_1"] <-"13"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "S1_2"] <-"14"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "S1_3"] <-"15"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "S2_1"] <-"16"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "S2_2"] <-"17"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "S2_3"] <-"18"

sample_df_G7$Tank[sample_df_G7$TreatmentTank == "C1_1"] <-"19"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "C1_2"] <-"20"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "C1_3"] <-"21"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "C2_1"] <-"22"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "C2_2"] <-"23"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "C2_3"] <-"24"

sample_df_G7$Tank[sample_df_G7$TreatmentTank == "M1_1"] <-"25"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "M1_2"] <-"26"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "M1_3"] <-"27"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "M2_1"] <-"28"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "M2_2"] <-"29"
sample_df_G7$Tank[sample_df_G7$TreatmentTank == "M2_3"] <-"30"

## Add CoralID column
CoralID <- paste0(sample_df_G7$TreatmentTank, "_", sample_df_G7$Fragment)
sample_df_G7 %>% dplyr::mutate(CoralID) -> sample_df_G7

## Add info for decontam and library size

### Add Sample_control column
sample_df_G7$Sample_control <- "True Sample"

### Add information to Sample_control column
sample_df_G7$Sample_control[sample_df_G7$Antibiotic == "Negative Control"] <- "Negative Control"

### Add DNA quantification column
readr::read_csv(here::here("Data/G7_QuantReading.csv")) -> quant_G7
sample_df_G7 <- merge(sample_df_G7, quant_G7, by = "Samples")

### Add LibrarySize columns (pre and post QC)
readr::read_csv(here::here("Data/00 - Read Preprocessing - Output/librarySize.csv")) -> librarySize_df
sample_df_G7 <- merge(sample_df_G7, librarySize_df, by = "Samples")

## Save metadata file 
saveRDS(sample_df_G7, here::here("Data/01 - Metadata - Output/sample_df_G7.rds"))
readr::write_csv(sample_df_G7, here::here("Data/01 - Metadata - Output/sample_df_G7.csv"))

## Supplemental Figure 2 (library size) ====
# Remove negative controls for plotting
subset(sample_df_G7, Antibiotic != "Negative Control") -> sample_df_G7
# Set order
sample_df_G7$Antibiotic <- factor(sample_df_G7$Antibiotic, c("Blank", "Ampicillin", "Streptomycin", "Ciprofloxacin", "Mixture"), ordered = TRUE)

# Pre QC
plot_preQC <-  ggplot(sample_df_G7, aes(x = as.factor(Time), y = PreQC, color = Antibiotic)) + 
  scale_color_manual(values = c("#666666", "#BC71BB", "#006DDB", "#008600", "#FF3D3D")) + geom_boxplot(lwd = 1.1, outlier.colour = "NA")

plot_preQC <- plot_preQC + geom_point(aes(color = Antibiotic), alpha = 0.5, 
                                          position = position_jitterdodge(jitter.width = 0.1)) +
  theme_bw() + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) 

plot_preQC <- plot_preQC + xlab("Time (Hours)") + ylab("Library Size") 
plot_preQC <- plot_preQC + labs(title = "Pre QC") 

# Post QC
plot_postQC <-  ggplot(sample_df_G7, aes(x = as.factor(Time), y = PostQC, color = Antibiotic)) + 
  scale_color_manual(values = c("#666666", "#BC71BB", "#006DDB", "#008600", "#FF3D3D")) + geom_boxplot(lwd = 1.1, outlier.colour = "NA")

plot_postQC <- plot_postQC + geom_point(aes(color = Antibiotic), alpha = 0.5, 
                                      position = position_jitterdodge(jitter.width = 0.1)) +
  theme_bw() + theme(axis.text = element_text(face = "bold", size = 11.5), 
                     axis.title = element_text(face = "bold", size = 12), title = element_text(face = "bold")) 

plot_postQC <- plot_postQC + xlab("Time (Hours)") + ylab("Library Size") 
plot_postQC <- plot_postQC + labs(title = "Post QC") 

librarySize <- ggarrange(plot_preQC, plot_postQC, nrow = 2, ncol = 1)

ggplot2::ggsave(here::here("Data/01 - Metadata - Output/librarySize.png"), librarySize,
               height = 350, width = 450, units = "mm",
              scale = 0.5, dpi = 1000)               

## Supplemental Figure 2 (library size) statistics ====

## Kruskal-Wallis
### Pre-QC; within-timepoint differences between antibiotic groups
broom::tidy(kruskal.test(PreQC ~ Antibiotic, data = subset(sample_df_G7, Time == "0"))) -> kruskal_PreQC_T0 # Not sig 
kruskal_PreQC_T0$Time <- "0"

broom::tidy(kruskal.test(PreQC ~ Antibiotic, data = subset(sample_df_G7, Time == "12"))) -> kruskal_PreQC_T12 # Sig
kruskal_PreQC_T12$Time <- "12"

broom::tidy(kruskal.test(PreQC ~ Antibiotic, data = subset(sample_df_G7, Time == "24"))) -> kruskal_PreQC_T24 # Sig 
kruskal_PreQC_T24$Time <- "24"

broom::tidy(kruskal.test(PreQC ~ Antibiotic, data = subset(sample_df_G7, Time == "48"))) -> kruskal_PreQC_T48 # Sig 
kruskal_PreQC_T48$Time <- "48"

broom::tidy(kruskal.test(PreQC ~ Antibiotic, data = subset(sample_df_G7, Time == "96"))) -> kruskal_PreQC_T96 # Sig 
kruskal_PreQC_T96$Time <- "96"

kruskal_PreQC <- bind_rows(kruskal_PreQC_T0, kruskal_PreQC_T12, kruskal_PreQC_T24, kruskal_PreQC_T48, kruskal_PreQC_T96)
kruskal_PreQC$QC <- "Pre-QC"

## Post hoc Dunn test 
dunn_Pre_T12 <- dunn_test(subset(sample_df_G7, Time == "12"), PreQC ~ Antibiotic, p.adjust.method = "BH")
dunn_Pre_T12$Time <- "12"

dunn_Pre_T24 <- dunn_test(subset(sample_df_G7, Time == "24"), PreQC ~ Antibiotic, p.adjust.method = "BH")
dunn_Pre_T24$Time <- "24"

dunn_Pre_T48 <- dunn_test(subset(sample_df_G7, Time == "48"), PreQC ~ Antibiotic, p.adjust.method = "BH")
dunn_Pre_T48$Time <- "48"

dunn_Pre_T96 <- dunn_test(subset(sample_df_G7, Time == "96"), PreQC ~ Antibiotic, p.adjust.method = "BH")
dunn_Pre_T96$Time <- "96"

dunn_preQC <- bind_rows(dunn_Pre_T12, dunn_Pre_T24, dunn_Pre_T48, dunn_Pre_T96)

## Post QC Kruskal-Wallis
broom::tidy(kruskal.test(PostQC ~ Antibiotic, data = subset(sample_df_G7, Time == "0"))) -> kruskal_PostQC_T0  # Not sig 
kruskal_PostQC_T0$Time <- "0"

broom::tidy(kruskal.test(PostQC ~ Antibiotic, data = subset(sample_df_G7, Time == "12"))) -> kruskal_PostQC_T12 # Sig 
kruskal_PostQC_T12$Time <- "12"

broom::tidy(kruskal.test(PostQC ~ Antibiotic, data = subset(sample_df_G7, Time == "24"))) -> kruskal_PostQC_T24 # Sig
kruskal_PostQC_T24$Time <- "24"

broom::tidy(kruskal.test(PostQC ~ Antibiotic, data = subset(sample_df_G7, Time == "48"))) -> kruskal_PostQC_T48 # Sig 
kruskal_PostQC_T48$Time <- "48"

broom::tidy(kruskal.test(PostQC ~ Antibiotic, data = subset(sample_df_G7, Time == "96"))) -> kruskal_PostQC_T96 # Sig
kruskal_PostQC_T96$Time <- "96"

kruskal_PostQC <- bind_rows(kruskal_PostQC_T0, kruskal_PostQC_T12, kruskal_PostQC_T24, kruskal_PostQC_T48, kruskal_PostQC_T96)
kruskal_PostQC$QC <- "Post-QC"

## Post hoc Dunn test 
dunn_Post_T12 <- dunn_test(subset(sample_df_G7, Time == "12"), PostQC ~ Antibiotic, p.adjust.method = "BH")
dunn_Post_T12$Time <- "12"

dunn_Post_T24 <- dunn_test(subset(sample_df_G7, Time == "24"), PostQC ~ Antibiotic, p.adjust.method = "BH")
dunn_Post_T24$Time <- "24"

dunn_Post_T48 <- dunn_test(subset(sample_df_G7, Time == "48"), PostQC ~ Antibiotic, p.adjust.method = "BH")
dunn_Post_T48$Time <- "48"

dunn_Post_T96 <- dunn_test(subset(sample_df_G7, Time == "96"), PostQC ~ Antibiotic, p.adjust.method = "BH")
dunn_Post_T96$Time <- "96"

dunn_postQC <- bind_rows(dunn_Post_T12, dunn_Post_T24, dunn_Post_T48, dunn_Post_T96)

## Save kruskal and dunn outputs
kruskal_librarySize <- bind_rows(kruskal_PreQC, kruskal_PostQC)
dunnTest_librarySize <- bind_rows(dunn_preQC, dunn_postQC)

## Save stats file 
write.csv(kruskal_librarySize, here::here("Data/01 - Metadata - Output/kruskal_libSize_stats.csv"))
write.csv(dunnTest_librarySize, here::here("Data/01 - Metadata - Output/dunnTest_libSize_stats.csv"))