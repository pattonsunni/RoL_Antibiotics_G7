# Author: Sunni Patton
# Last edited: 09/16/2024
# Title: Network Analysis - SpiecEasi MicroEco
# Adapted from: https://chiliubio.github.io/microeco_tutorial/basic-class.html

## Set seed ====
set.seed(123)

## Load libraries ====
library(phyloseq)
library(here)
library(microViz)
library(file2meco) # used for converting a phyloseq object to microtable
library(magrittr)
library(meconetcomp)
library(microeco)
library(SpiecEasi)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(ggpubr)
library(stringr)
library(igraph)
library(MicEco) # for venn diagrams
library(seqateurs)
library(RCy3)
library(httr)
library(patchwork)

## Load data ====
readRDS(here::here("Data/05 - Phyloseq - Output/ps.Tree.rds")) -> ps

## Fix data 
tax_fix(
  ps,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
) -> ps

# Rename taxa
taxa_change <- data.frame(tax_table(ps))

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
tax_table(ps) <- as.matrix(taxa_change)

# Add additional column 
tax.table <- tax_table(ps)
as.data.frame(tax.table) -> tax.table
tax.table$ASV2 <- paste0(tax.table$Genus, " (", tax.table$ASV, ")")
tax.table <- tax_table(as.matrix(tax.table))
merge_phyloseq(ps, tax.table) -> ps
ps@tax_table <- ps@tax_table[,c(1:7, 9, 8)]
colnames(ps@tax_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV", "ASV2")

## Subset data ====
subset_samples(ps, Time != "0") -> ps

# Remove sample G7_T96_S1_3_d (contamination) --> also need to fix the total samples figure
subset_samples(ps, Samples != "G7_T96_S1_3_d") -> ps

ps.blank <- subset_samples(ps, Antibiotic == "Blank")
ps.blank <- prune_taxa(taxa_sums(ps.blank@otu_table) > 0, ps.blank) 

ps.mix <- subset_samples(ps, Antibiotic == "Mixture")
ps.mix <- prune_taxa(taxa_sums(ps.mix@otu_table) > 0, ps.mix) 

ps.amp <- subset_samples(ps, Antibiotic == "Ampicillin")
ps.amp <- prune_taxa(taxa_sums(ps.amp@otu_table) > 0, ps.amp) 

ps.strep <- subset_samples(ps, Antibiotic == "Streptomycin")
ps.strep <- prune_taxa(taxa_sums(ps.strep@otu_table) > 0, ps.strep) 

ps.cip <- subset_samples(ps, Antibiotic == "Ciprofloxacin")
ps.cip <- prune_taxa(taxa_sums(ps.cip@otu_table) > 0, ps.cip) 

## Blank ====
pst <- fast_melt(ps.blank) # Takes a phyloseq object and returns a melted dataframe with sample ID, variable name, and abundance value for each observation

prevdt <- pst[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)),
              by = taxaID]

keepTaxa <- prevdt[(Prevalence > 1), taxaID] # replace with appropriate number for each network
ps.blank <- prune_taxa(keepTaxa, ps.blank) 
ntaxa(ps.blank)

## check
pst <- fast_melt(ps.blank) 

prevdt <- pst[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)),
              by = taxaID]

# Blank: 191 taxa

meco_blank <- phyloseq2meco(ps.blank) # creates microtable R6 object

# Make copy of data
tmp <- clone(meco_blank)

# Trim all files in the object
tmp$tidy_dataset() 

# Network construction 
t1_blank <- trans_network$new(dataset = tmp, cor_method = NULL, taxa_level = "ASV")
## Blank: 191 features remain

# Calculate network (Creates 'res_network' which is igraph object)
t1_blank$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", lambda.min.ratio=1e-2, pulsar.params=list(rep.num=50, seed = 10010)) # utilizing SpiecEasi method with meinshausen-buhlmann's neighborhood selection as the estimation method

# Get edge table
t1_blank$get_edge_table()

# Calculate modules
t1_blank$cal_module(method = "cluster_fast_greedy") # identifies 23 modules
## Blank: 23 modules 

# Get node table
t1_blank$get_node_table()

# Get adjacency matrix
t1_blank$get_adjacency_matrix()

# Calculate network attributes
t1_blank$cal_network_attr() # saves network properties

## Calculate edge sums
t1_blank$cal_sum_links(taxa_level = "ASV") # Sums the total number of connections (edges) for each taxon at the level specified. This is then used in a chord diagram or alluvial plot

# Save network 
## gexf for visualization in Gephi
t1_blank$save_network(filepath = here::here("Data/10 - Network Analysis - Output/blank_gephi_new.gexf"))
saveRDS(t1_blank, here::here("Data/10 - Network Analysis - Output/t1_blank_test_new.rds"))

## Visualization in cytoscape
createNetworkFromIgraph(t1_blank$res_network,"blank_network_cyto_new") # must have cytoscape running

## Save node table
write.csv(t1_blank$res_node_table, here::here("Data/10 - Network Analysis - Output/blank_nodeTable.csv"))

## Prepare data for alluvial plot
### Save cal_sum_links output
t1_blank$res_sum_links_pos -> blank_pos
as.table(blank_pos) -> blank_pos
as.data.frame(blank_pos) -> blank_pos
blank_pos$Direction <- "Positive"

t1_blank$res_sum_links_neg -> blank_neg
as.table(blank_neg) -> blank_neg
as.data.frame(blank_neg) -> blank_neg
blank_neg$Direction <- "Negative"

## Remove rows that don't contain a link (i.e., freq = 0)
subset(blank_pos, Freq > 0) -> blank_pos
subset(blank_neg, Freq > 0) -> blank_neg


subset(blank_pos, blank_pos$Var1 == "a__Aquarickettsia (ASV6)" |
         blank_pos$Var1 == "a__Campylobacterales Order (ASV1)" |
         blank_pos$Var1 == "a__Alteromonas (ASV2)" |
         blank_pos$Var1 == "a__Helicobacteraceae Family (ASV3)") -> blank_pos_subset

subset(blank_neg, blank_neg$Var1 == "a__Aquarickettsia (ASV6)" |
         blank_neg$Var1 == "a__Campylobacterales Order (ASV1)" |
         blank_neg$Var1 == "a__Alteromonas (ASV2)" |
         blank_neg$Var1 == "a__Helicobacteraceae Family (ASV3)") -> blank_neg_subset
### Clean up phyla names
blank_pos_subset$Var1 <- str_remove_all(blank_pos_subset$Var1, "a__")
blank_pos_subset$Var2 <- str_remove_all(blank_pos_subset$Var2, "a__")
blank_neg_subset$Var1 <- str_remove_all(blank_neg_subset$Var1, "a__")
blank_neg_subset$Var2 <- str_remove_all(blank_neg_subset$Var2, "a__")


colors <- c("#A6CEE2", "#1E78B4", "#B2DF8A", "#33A12C", "#FB9999", "#E3201D", "#FDBF6F", "#FF7E02", "#CAB2D6", "#6A3D9A", "#E7AB02", "#B15929", "#878787", "#4D4D4D", "#C5247E", "#8F1652", "#FFFD9A", "cyan")

pos_blank <- ggplot(data = blank_pos_subset,
                    aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) + geom_stratum() +
  geom_text(stat = "stratum", size=2.6, fontface = "bold",
            aes(label = after_stat(stratum))) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values = colors) + facet_wrap(~Direction) +
  scale_x_discrete(limits = c("Source", "Target"),
                   expand = c(0.15, 0.05)) + theme(axis.title.y = element_blank()) + 
  theme(axis.text = element_text(face = "bold", size = 10.5), axis.title = element_text(face = "bold", size = 12)) +
  theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_rect(fill = "#24B751")) + ggtitle("") + theme(plot.title = element_text(face="bold", color="#008600", size = 12)) + scale_y_continuous(breaks = seq(0, 8, by = 2))


neg_blank <- ggplot(data = blank_neg_subset,
                    aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) + geom_stratum() +
  geom_text(stat = "stratum", size=2.6, fontface = "bold",
            aes(label = after_stat(stratum), )) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values = colors) + facet_wrap(~Direction) +
  scale_x_discrete(limits = c("Source", "Target"),
                   expand = c(0.15, 0.05)) + ylab("Frequency") + 
  theme(axis.text = element_text(face = "bold", size = 10.5), axis.title = element_text(face = "bold", size = 12)) +
  theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_rect(fill = "#FA6F6F")) + ggtitle("Blank") + theme(plot.title = element_text(face="bold", color="black", size = 14)) + scale_y_continuous(breaks = seq(0, 2, by = 1))

ggarrange(neg_blank, pos_blank) -> plot_blank

## Calculate metrics 
View(t1_blank$res_node_table) # already calculated some metrics in node table


View(prevdt)

t1_blank$trans_comm(use_col = "module", abundance = FALSE)


## Mixture ====
pst <- fast_melt(ps.mix) # Takes a phyloseq object and returns a melted dataframe with sample ID, variable name, and abundance value for each observation

prevdt <- pst[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)),
              by = taxaID]

keepTaxa <- prevdt[(Prevalence > 1), taxaID] # replace with appropriate number for each network
ps.mix <- prune_taxa(keepTaxa, ps.mix) 
ntaxa(ps.mix)

# Blank: 191 taxa
# Mix: 127 taxa

meco_mix <- phyloseq2meco(ps.mix) # creates microtable R6 object

# Make copy of data
tmp <- clone(meco_mix)

# Trim all files in the object
tmp$tidy_dataset() 

# Network construction 
t1_mix <- trans_network$new(dataset = tmp, cor_method = NULL, taxa_level = "ASV")
## Blank: 191 features remain

# Calculate network (Creates 'res_network' which is igraph object)
t1_mix$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", lambda.min.ratio=1e-2, pulsar.params=list(rep.num=50, seed = 10010)) # utilizing SpiecEasi method with meinshausen-buhlmann's neighborhood selection as the estimation method

# Get edge table
t1_mix$get_edge_table()

# Calculate modules
t1_mix$cal_module(method = "cluster_fast_greedy") 
## Blank: 23 modules 
## Mix: 17 modules

# Get node table
t1_mix$get_node_table()

# Get adjacency matrix
t1_mix$get_adjacency_matrix()

# Calculate network attributes
t1_mix$cal_network_attr() # saves network properties

## Calculate edge sums
t1_mix$cal_sum_links(taxa_level = "ASV") # Sums the total number of connections (edges) for each taxon at the level specified. This is then used in a chord diagram or alluvial plot
# no negative edges found

# Save network 
## gexf for visualization in Gephi
t1_mix$save_network(filepath = here::here("Data/10 - Network Analysis - Output/mix_gephi_new.gexf"))
saveRDS(t1_mix, here::here("Data/10 - Network Analysis - Output/t1_mix_new.rds"))

## Visualization in cytoscape
createNetworkFromIgraph(t1_mix$res_network,"mix_network_cyto_new") # must have cytoscape running

## Save node table
write.csv(t1_mix$res_node_table, here::here("Data/10 - Network Analysis - Output/mix_nodeTable.csv"))

## Prepare data for alluvial plot
### Save cal_sum_links output
t1_mix$res_sum_links_pos -> mix_pos
as.table(mix_pos) -> mix_pos
as.data.frame(mix_pos) -> mix_pos
mix_pos$Direction <- "Positive"

## Remove rows that don't contain a link (i.e., freq = 0)
subset(mix_pos, Freq > 0) -> mix_pos

subset(mix_pos, mix_pos$Var1 == "a__Campylobacterales Order (ASV1)" |
         mix_pos$Var1 == "a__Alteromonas (ASV2)" |
         mix_pos$Var1 == "a__Helicobacteraceae Family (ASV3)") -> mix_pos_subset

### Clean up phyla names
mix_pos_subset$Var1 <- str_remove_all(mix_pos_subset$Var1, "a__")
mix_pos_subset$Var2 <- str_remove_all(mix_pos_subset$Var2, "a__")

## Add species name if known
mix_pos_subset$Var2[mix_pos_subset$Var2 == "Idiomarina (ASV38)"] <- "Idiomarina andamanensis (ASV38)"

colors <- c("#A6CEE2", "#1E78B4", "#B2DF8A", "#33A12C", "#FB9999", "#E3201D", "#FDBF6F", "#FF7E02", "#CAB2D6", "#6A3D9A", "#E7AB02", "#B15929", "#878787", "#4D4D4D", "#C5247E", "#8F1652", "#FFFD9A", "cyan")

pos_mix <- ggplot(data = mix_pos_subset,
                  aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) + geom_stratum() +
  geom_text(stat = "stratum", size=2.6, fontface = "bold",
            aes(label = after_stat(stratum))) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values = colors) + facet_wrap(~Direction) +
  scale_x_discrete(limits = c("Source", "Target"),
                   expand = c(0.15, 0.05)) + ylab("Frequency") +
  theme(axis.text = element_text(face = "bold", size = 10.5), axis.title = element_text(face = "bold", size = 12)) +
  theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_rect(fill = "#24B751")) + ggtitle("Mixture") + theme(plot.title = element_text(face="bold", color="#FF3D3D", size = 14)) #+ scale_y_continuous(breaks = seq(0, 10, by = 2))

ggarrange(pos_mix) -> plot_mix

## Ampicillin ====
pst <- fast_melt(ps.amp) # Takes a phyloseq object and returns a melted dataframe with sample ID, variable name, and abundance value for each observation

prevdt <- pst[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)),
              by = taxaID]

keepTaxa <- prevdt[(Prevalence > 1), taxaID] # replace with appropriate number for each network
ps.amp <- prune_taxa(keepTaxa, ps.amp) 
ntaxa(ps.amp)

# Blank: 191 taxa
# Mix: 127 taxa
# Amp: 197 taxa

meco_amp <- phyloseq2meco(ps.amp) # creates microtable R6 object

# Make copy of data
tmp <- clone(meco_amp)

# Trim all files in the object
tmp$tidy_dataset() 

# Network construction 
t1_amp <- trans_network$new(dataset = tmp, cor_method = NULL, taxa_level = "ASV")
## Blank: 191 features remain

# Calculate network (Creates 'res_network' which is igraph object)
t1_amp$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", lambda.min.ratio=1e-2, pulsar.params=list(rep.num=50, seed = 10010)) # utilizing SpiecEasi method with meinshausen-buhlmann's neighborhood selection as the estimation method

# Get edge table
t1_amp$get_edge_table()

# Calculate modules
t1_amp$cal_module(method = "cluster_fast_greedy") 
## Blank: 23 modules 
## Amp: 21 modules

# Get node table
t1_amp$get_node_table()

# Get adjacency matrix
t1_amp$get_adjacency_matrix()

# Calculate network attributes
t1_amp$cal_network_attr() # saves network properties

## Calculate edge sums
t1_amp$cal_sum_links(taxa_level = "ASV") # Sums the total number of connections (edges) for each taxon at the level specified. This is then used in a chord diagram or alluvial plot

# Save network 
## gexf for visualization in Gephi
t1_amp$save_network(filepath = here::here("Data/10 - Network Analysis - Output/amp_gephi_new.gexf"))
saveRDS(t1_amp, here::here("Data/10 - Network Analysis - Output/t1_amp_new.rds"))

## Visualization in cytoscape
createNetworkFromIgraph(t1_amp$res_network,"amp_network_cyto_new") # must have cytoscape running

## Save node table
write.csv(t1_amp$res_node_table, here::here("Data/10 - Network Analysis - Output/amp_nodeTable.csv"))

## Prepare data for alluvial plot
### Save cal_sum_links output
t1_amp$res_sum_links_pos -> amp_pos
as.table(amp_pos) -> amp_pos
as.data.frame(amp_pos) -> amp_pos
amp_pos$Direction <- "Positive"

t1_amp$res_sum_links_neg -> amp_neg
as.table(amp_neg) -> amp_neg
as.data.frame(amp_neg) -> amp_neg
amp_neg$Direction <- "Negative"

## Remove rows that don't contain a link (i.e., freq = 0)
subset(amp_pos, Freq > 0) -> amp_pos
subset(amp_neg, Freq > 0) -> amp_neg

subset(amp_pos, amp_pos$Var1 == "a__Campylobacterales Order (ASV1)" |
         amp_pos$Var1 == "a__Alteromonas (ASV2)" |
         amp_pos$Var1 == "a__Helicobacteraceae Family (ASV3)") -> amp_pos_subset

subset(amp_neg, amp_neg$Var1 == "a__Campylobacterales Order (ASV1)" |
         amp_neg$Var1 == "a__Alteromonas (ASV2)" |
         amp_neg$Var1 == "a__Helicobacteraceae Family (ASV3)") -> amp_neg_subset

### Clean up phyla names
amp_pos_subset$Var1 <- str_remove_all(amp_pos_subset$Var1, "a__")
amp_pos_subset$Var2 <- str_remove_all(amp_pos_subset$Var2, "a__")
amp_neg_subset$Var1 <- str_remove_all(amp_neg_subset$Var1, "a__")
amp_neg_subset$Var2 <- str_remove_all(amp_neg_subset$Var2, "a__")

colors <- c("#A6CEE2", "#1E78B4", "#B2DF8A", "#33A12C", "#FB9999", "#E3201D", "#FDBF6F", "#FF7E02", "#CAB2D6", "#6A3D9A", "#E7AB02", "#B15929", "#878787", "#4D4D4D", "#C5247E", "#8F1652", "#FFFD9A", "cyan")

pos_amp <- ggplot(data = amp_pos_subset,
                  aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) + geom_stratum() +
  geom_text(stat = "stratum", size=2.6, fontface = "bold",
            aes(label = after_stat(stratum))) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values = colors) + facet_wrap(~Direction) +
  scale_x_discrete(limits = c("Source", "Target"),
                   expand = c(0.15, 0.05)) + theme(axis.title.y = element_blank()) + 
  theme(axis.text = element_text(face = "bold", size = 10.5), axis.title = element_text(face = "bold", size = 12)) +
  theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_rect(fill = "#24B751")) + ggtitle("") + theme(plot.title = element_text(face="bold", color="#008600", size = 12))


neg_amp <- ggplot(data = amp_neg_subset,
                  aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) + geom_stratum() +
  geom_text(stat = "stratum", size=2.6, fontface = "bold",
            aes(label = after_stat(stratum), )) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values = colors) + facet_wrap(~Direction) +
  scale_x_discrete(limits = c("Source", "Target"),
                   expand = c(0.15, 0.05)) + ylab("Frequency") + 
  theme(axis.text = element_text(face = "bold", size = 10.5), axis.title = element_text(face = "bold", size = 12)) +
  theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_rect(fill = "#FA6F6F")) + ggtitle("Ampicillin") + theme(plot.title = element_text(face="bold", color="#BC71BB", size = 14)) 

ggarrange(neg_amp, pos_amp) -> plot_amp

## Streptomycin ====
pst <- fast_melt(ps.strep) # Takes a phyloseq object and returns a melted dataframe with sample ID, variable name, and abundance value for each observation

prevdt <- pst[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)),
              by = taxaID]

keepTaxa <- prevdt[(Prevalence > 1), taxaID] # replace with appropriate number for each network
ps.strep <- prune_taxa(keepTaxa, ps.strep) 
ntaxa(ps.strep)

# Blank: 191 taxa
# Mix: 127 taxa
# Amp: 197 taxa
# Strep: 188 taxa

meco_strep <- phyloseq2meco(ps.strep) # creates microtable R6 object

# Make copy of data
tmp <- clone(meco_strep)

# Trim all files in the object
tmp$tidy_dataset() 

# Network construction 
t1_strep <- trans_network$new(dataset = tmp, cor_method = NULL, taxa_level = "ASV")
## Blank: 191 features remain

# Calculate network (Creates 'res_network' which is igraph object)
t1_strep$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", lambda.min.ratio=1e-2, pulsar.params=list(rep.num=50, seed = 10010)) # utilizing SpiecEasi method with meinshausen-buhlmann's neighborhood selection as the estimation method

# Get edge table
t1_strep$get_edge_table()

# Calculate modules
t1_strep$cal_module(method = "cluster_fast_greedy") 
## Blank: 23 modules 
## Amp: 21 modules
## Strep: 20 modules

# Get node table
t1_strep$get_node_table()

# Get adjacency matrix
t1_strep$get_adjacency_matrix()

# Calculate network attributes
t1_strep$cal_network_attr() # saves network properties

## Calculate edge sums
t1_strep$cal_sum_links(taxa_level = "ASV") # Sums the total number of connections (edges) for each taxon at the level specified. This is then used in a chord diagram or alluvial plot

# Save network 
## gexf for visualization in Gephi
t1_strep$save_network(filepath = here::here("Data/10 - Network Analysis - Output/strep_gephi_new.gexf"))
saveRDS(t1_strep, here::here("Data/10 - Network Analysis - Output/t1_strep_new.rds"))

## Visualization in cytoscape
createNetworkFromIgraph(t1_strep$res_network,"strep_network_cyto_new") # must have cytoscape running

## Save node table
write.csv(t1_strep$res_node_table, here::here("Data/10 - Network Analysis - Output/strep_nodeTable.csv"))

## Prepare data for alluvial plot
### Save cal_sum_links output
t1_strep$res_sum_links_pos -> strep_pos
as.table(strep_pos) -> strep_pos
as.data.frame(strep_pos) -> strep_pos
strep_pos$Direction <- "Positive"

t1_strep$res_sum_links_neg -> strep_neg
as.table(strep_neg) -> strep_neg
as.data.frame(strep_neg) -> strep_neg
strep_neg$Direction <- "Negative"

## Remove rows that don't contain a link (i.e., freq = 0)
subset(strep_pos, Freq > 0) -> strep_pos
subset(strep_neg, Freq > 0) -> strep_neg

subset(strep_pos, strep_pos$Var1 == "a__Campylobacterales Order (ASV1)" |
         strep_pos$Var1 == "a__Alteromonas (ASV2)" |
         strep_pos$Var1 == "a__Helicobacteraceae Family (ASV3)") -> strep_pos_subset

subset(strep_neg, strep_neg$Var1 == "a__Campylobacterales Order (ASV1)" |
         strep_neg$Var1 == "a__Alteromonas (ASV2)" |
         strep_neg$Var1 == "a__Helicobacteraceae Family (ASV3)") -> strep_neg_subset

### Clean up phyla names
strep_pos_subset$Var1 <- str_remove_all(strep_pos_subset$Var1, "a__")
strep_pos_subset$Var2 <- str_remove_all(strep_pos_subset$Var2, "a__")
strep_neg_subset$Var1 <- str_remove_all(strep_neg_subset$Var1, "a__")
strep_neg_subset$Var2 <- str_remove_all(strep_neg_subset$Var2, "a__")

## Add species name if known
strep_neg_subset$Var2[strep_neg_subset$Var2 == "Winogradskyella (ASV25)"] <- "Winogradskyella poriferorum (ASV25)"

colors <- c("#A6CEE2", "#1E78B4", "#B2DF8A", "#33A12C", "#FB9999", "#E3201D", "#FDBF6F", "#FF7E02", "#CAB2D6", "#6A3D9A", "#E7AB02", "#B15929", "#878787", "#4D4D4D", "#C5247E", "#8F1652", "#FFFD9A", "cyan")

pos_strep <- ggplot(data = strep_pos_subset,
                    aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) + geom_stratum() +
  geom_text(stat = "stratum", size=2.6, fontface = "bold",
            aes(label = after_stat(stratum))) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values = colors) + facet_wrap(~Direction) +
  scale_x_discrete(limits = c("Source", "Target"),
                   expand = c(0.15, 0.05)) + theme(axis.title.y = element_blank()) + 
  theme(axis.text = element_text(face = "bold", size = 10.5), axis.title = element_text(face = "bold", size = 12)) +
  theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_rect(fill = "#24B751")) + ggtitle("Streptomycin") + theme(plot.title = element_text(face="bold", color="#006DDB", size = 12)) + scale_y_continuous(breaks = seq(0, 10, by = 2))

ggarrange(pos_strep) -> plot_strep

## Ciprofloxacin ====
pst <- fast_melt(ps.cip) # Takes a phyloseq object and returns a melted dataframe with sample ID, variable name, and abundance value for each observation

prevdt <- pst[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)),
              by = taxaID]

keepTaxa <- prevdt[(Prevalence > 1), taxaID] # replace with appropriate number for each network
ps.cip <- prune_taxa(keepTaxa, ps.cip) 
ntaxa(ps.cip)

# Blank: 191 taxa
# Mix: 127 taxa
# Amp: 197 taxa
# Strep: 188 taxa
# Cip: 195 taxa

meco_cip <- phyloseq2meco(ps.cip) # creates microtable R6 object

# Make copy of data
tmp <- clone(meco_cip)

# Trim all files in the object
tmp$tidy_dataset() 

# Network construction 
t1_cip <- trans_network$new(dataset = tmp, cor_method = NULL, taxa_level = "ASV")
## Blank: 191 features remain

# Calculate network (Creates 'res_network' which is igraph object)
t1_cip$cal_network(network_method = "SpiecEasi", SpiecEasi_method = "mb", lambda.min.ratio=1e-2, pulsar.params=list(rep.num=50, seed = 10010)) # utilizing SpiecEasi method with meinshausen-buhlmann's neighborhood selection as the estimation method

# Get edge table
t1_cip$get_edge_table()

# Calculate modules
t1_cip$cal_module(method = "cluster_fast_greedy") 
## Blank: 23 modules 
## Mix: 17 modules 
## Amp: 21 modules
## Strep: 20 modules
## Cip: 16 modules

# Get node table
t1_cip$get_node_table()

# Get adjacency matrix
t1_cip$get_adjacency_matrix()

# Calculate network attributes
t1_cip$cal_network_attr() # saves network properties

## Calculate edge sums
t1_cip$cal_sum_links(taxa_level = "ASV") # Sums the total number of connections (edges) for each taxon at the level specified. This is then used in a chord diagram or alluvial plot

# Save network 
## gexf for visualization in Gephi
t1_cip$save_network(filepath = here::here("Data/10 - Network Analysis - Output/cip_gephi_new.gexf"))
saveRDS(t1_cip, here::here("Data/10 - Network Analysis - Output/t1_cip_new.rds"))

## Visualization in cytoscape
createNetworkFromIgraph(t1_cip$res_network,"cip_network_cyto_new") # must have cytoscape running

## Save node table 
write.csv(t1_cip$res_node_table, here::here("Data/10 - Network Analysis - Output/cip_nodeTable.csv"))

## Prepare data for alluvial plot
### Save cal_sum_links output
t1_cip$res_sum_links_pos -> cip_pos
as.table(cip_pos) -> cip_pos
as.data.frame(cip_pos) -> cip_pos
cip_pos$Direction <- "Positive"

t1_cip$res_sum_links_neg -> cip_neg
as.table(cip_neg) -> cip_neg
as.data.frame(cip_neg) -> cip_neg
cip_neg$Direction <- "Negative"

## Remove rows that don't contain a link (i.e., freq = 0)
subset(cip_pos, Freq > 0) -> cip_pos
subset(cip_neg, Freq > 0) -> cip_neg

subset(cip_pos, cip_pos$Var1 == "a__Campylobacterales Order (ASV1)" |
         cip_pos$Var1 == "a__Alteromonas (ASV2)" |
         cip_pos$Var1 == "a__Helicobacteraceae Family (ASV3)") -> cip_pos_subset

subset(cip_neg, cip_neg$Var1 == "a__Campylobacterales Order (ASV1)" |
         cip_neg$Var1 == "a__Alteromonas (ASV2)" |
         cip_neg$Var1 == "a__Helicobacteraceae Family (ASV3)") -> cip_neg_subset

### Clean up phyla names
cip_pos_subset$Var1 <- str_remove_all(cip_pos_subset$Var1, "a__")
cip_pos_subset$Var2 <- str_remove_all(cip_pos_subset$Var2, "a__")
cip_neg_subset$Var1 <- str_remove_all(cip_neg_subset$Var1, "a__")
cip_neg_subset$Var2 <- str_remove_all(cip_neg_subset$Var2, "a__")

colors <- c("#A6CEE2", "#1E78B4", "#B2DF8A", "#33A12C", "#FB9999", "#E3201D", "#FDBF6F", "#FF7E02", "#CAB2D6", "#6A3D9A", "#E7AB02", "#B15929", "#878787", "#4D4D4D", "#C5247E", "#8F1652", "#FFFD9A", "cyan")

pos_cip <- ggplot(data = cip_pos_subset,
                  aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) + geom_stratum() +
  geom_text(stat = "stratum", size=2.6, fontface = "bold",
            aes(label = after_stat(stratum))) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values = colors) + facet_wrap(~Direction) +
  scale_x_discrete(limits = c("Source", "Target"),
                   expand = c(0.15, 0.05)) + theme(axis.title.y = element_blank()) + 
  theme(axis.text = element_text(face = "bold", size = 10.5), axis.title = element_text(face = "bold", size = 12)) +
  theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_rect(fill = "#24B751")) + ggtitle("") + theme(plot.title = element_text(face="bold", color="#008600", size = 12)) + scale_y_continuous(breaks = seq(0, 2, by = 1))


neg_cip <- ggplot(data = cip_neg_subset,
                  aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) + geom_stratum() +
  geom_text(stat = "stratum", size=2.6, fontface = "bold",
            aes(label = after_stat(stratum), )) + theme_bw() + theme(legend.position = "none") + 
  scale_fill_manual(values = colors) + facet_wrap(~Direction) +
  scale_x_discrete(limits = c("Source", "Target"),
                   expand = c(0.15, 0.05)) + ylab("Frequency") + 
  theme(axis.text = element_text(face = "bold", size = 10.5), axis.title = element_text(face = "bold", size = 12)) +
  theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_rect(fill = "#FA6F6F")) + ggtitle("Ciprofloxacin") + theme(plot.title = element_text(face="bold", color="#008600", size = 14)) + scale_y_continuous(breaks = seq(0, 2, by = 1))

ggarrange(neg_cip, pos_cip) -> plot_cip

## Arrange ====

plot_blank / (plot_mix | plot_strep) / plot_amp / plot_cip -> plot_all

ggarrange(plot_blank, plot_mix, plot_amp, plot_strep, plot_cip, nrow =5) -> plot_all


ggplot2::ggsave(here::here("Data/10 - Network Analysis - Output/alluvial_ASV_all.png"), plot_all,
                height = 750, width = 900, units = "mm",
                scale = 0.5, dpi = 1000)



## Network properties
### Clustering coefficient/transitivity
transitivity(t1_blank$res_network, type = "global") # type = c("weighted", "global")
### Blank: 0.1757, Mix: 0.1355, Amp: 0.17, Strep:0.215, Cip: 0.1889

### Modularity
wt_cip <- walktrap.community(t1_cip$res_network)
modularity(t1_cip$res_network, membership(wt_cip))

### Blank: 0.54, Mix: 0.558, Amp: 0.525, Strep: 0.43, Cip: 0.46

