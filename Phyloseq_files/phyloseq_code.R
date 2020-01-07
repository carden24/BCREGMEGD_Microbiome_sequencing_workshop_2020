# I recommend to run this to remove previous results just in case
rm(list=ls()) 


# This is my working directory, change it to something more useful for you

setwd("E:/Libraries/Dropbox/tutorial/") #Windows examples
setwd("~/pseq")       # Linux example
getwd() # This command reports the current working directory

library(phyloseq)
library(ggplot2) # This package is required for data visualization
library(vegan) # This package is required for the diversity analysis
library(tidyr) # These two package help with data manipulation
library(dplyr)
source("miseqR.R")


# Create variables for the files we will import from Mothur
shared_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared"
tax_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy"
metadata_file = "mouse.dpw.metadata"
biom_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.biom"
tree_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.phylip.tre"

# Import from the biom file
data = import_biom(biom_file)

data
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 361 taxa and 8 samples ]
# sample_data() Sample Data:       [ 8 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 361 taxa by 6 taxonomic ranks ]



# Importing Mothur files
mothur_data0 <- import_mothur(mothur_shared_file = shared_file,
                              mothur_constaxonomy_file = tax_file)

# Adding metadata

metadata = read.delim(metadata_file, header=T, row.names=1)
metadata = sample_data(metadata) # Converts the metadata table into a metadata table ready for *Phyloseq*
mothur_data = merge_phyloseq(mothur_data0, metadata)

mothur_data
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 361 taxa and 8 samples ]
# sample_data() Sample Data:       [ 8 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 361 taxa by 6 taxonomic ranks ]


### Accessing components

mothur_data
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 361 taxa and 8 samples ]
# sample_data() Sample Data:       [ 8 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 361 taxa by 6 taxonomic ranks ]

# Basic functions

otu_table(mothur_data)		# Reports the OTU table
sample_data(mothur_data)	# Reports the sample information such as experimental design or metadata
tax_table(mothur_data)		# Reports the taxonomic classification of each OTUs

ntaxa(mothur_data) 			# Number of OTUs
# [1] 361

nsamples(mothur_data) 		# Number of samples
# [1] 8

sample_names(mothur_data)	# Sample names
# [1] "F3D0"   "F3D1"   "F3D141" "F3D142" "F3D143" "F3D144" "F3D2"   "F3D3"  

taxa_names(mothur_data)		# Taxa names

taxa_sums(mothur_data)		# Sum of abundances or an OTU for all the samples


# Replacing colnames of taxa 
colnames(tax_table(mothur_data)) # Reports the names of the taxonomic ranks
#[1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6"

# Replacing the names
colnames(tax_table(mothur_data)) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family", "Genus")

colnames(tax_table(mothur_data))
# [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus" 





data("GlobalPatterns")
GlobalPatterns

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
# sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
# tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]

data("soilrep")
soilrep

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 16825 taxa and 56 samples ]
# sample_data() Sample Data:       [ 56 samples by 4 sample variables ]


# Sample distribution

sample_sum_df <- data.frame(sum = sample_sums(soilrep)) 

# Exploring mean, maximum, and minimum
min(sample_sum_df$sum)
# [1] 889

median(sample_sum_df$sum)
# [1] 1558.5

max(sample_sum_df$sum)
# [1] 4352

# Reads per samples histogram
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color="Black", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts")  + ylab("Frequency")


# Selecting and filtering samples


## To remove samples with fewer than 300 reads
mothur_data2 = prune_samples(names(which(sample_sums(mothur_data) >= 3000)), mothur_data)

mothur_data # <- Original dataset
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 522 taxa and 19 samples ]
# sample_data() Sample Data:       [ 19 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 522 taxa by 6 taxonomic ranks ]

mothur_data2 # <- Filtered dataset
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 522 taxa and 17 samples ]
# sample_data() Sample Data:       [ 17 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 522 taxa by 6 taxonomic ranks ]


# Selecting samples according to their metadata, in this case we want only

# samples that that "soil" as the "SampleType"
# In this example we further filter OTUs that have no data after removing those
# samples with prune_taxa() (see below)
global_soil <- GlobalPatterns %>%
  subset_samples(SampleType == "Soil") %>%
  prune_taxa(taxa_sums(.) > 0, .)




# Selecting and filtering OTUs

# The following line removes OTUs that have five or fewer reads
# in the whole dataset
mothur_data3 = prune_taxa(names(which(taxa_sums(mothur_data2) >= 5)), mothur_data2)

mothur_data3
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 174 taxa and 6 samples ]
# sample_data() Sample Data:       [ 6 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 174 taxa by 6 taxonomic ranks ]



# Subsetting taxa

# This will select all OTUs whose Phylum is Actinobacteria
Actinobacteria = subset_taxa(GlobalPatterns, Phylum == "Actinobacteria")

# This will get all Pseudomonas by selecting OTUs whose Genus is Pseudomonas
Pseudomonas = subset_taxa(GlobalPatterns, Genus == "Pseudomonas")



# Normalization

# This will take a phyloseq object and transform counts into relative abundance (0 - 1 range)
rel_abun_soil =  transform_sample_counts(soilrep, function(x) {x/sum(x)}) 

# This will do a square-root transformation of the counts
rel_abun_sqrt_soil =  transform_sample_counts(rel_abun_soil, function(x) sqrt(x))

set.seed(10) # Set the random seed to make the work reproducible
rare_soil <-rarefy_even_depth(soilrep, rngseed=TRUE)



# Community visualization


## Creating a phyla level plot

### Collapse data at Phylum level
global_phyla = tax_glom(physeq = GlobalPatterns, taxrank = "Phylum")

### Converting to relative data
global_phyla_rel_abu = transform_sample_counts(physeq = global_phyla, function(x){x/sum(x)})
plot_bar(global_phyla_rel_abu, fill="Phylum")


# This removes groups with abundances smaller than 5% 
global_phyla_rel_abu_clean = prune_taxa(taxa_sums(global_phyla_rel_abu) > 0.05, global_phyla_rel_abu)

phyla_plot = plot_bar(global_phyla_rel_abu_clean, fill="Phylum")
phyla_plot 
phyla_plot + facet_wrap(~SampleType, nrow=1, scales="free_x")


Acidos = subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
Acidos2 = tax_glom(Acidos, taxrank = "Class")
Acidos_rel_abu = transform_sample_counts(Acidos2, function(x){x/sum(x)})
Acido_plot = plot_bar(Acidos_rel_abu, fill="Class")
Acido_plot



# Heatmaps

# A heatmap without sample reordering
plot_heatmap(Acidos_rel_abu,  method=NULL, sample.label="SampleType", sample.order="SampleType")


# A heatmap with the samples reordered according to their sample similarity
plot_heatmap(Acidos_rel_abu,  # Object to plot
             method="NMDS",   # Method for the ordination
             distance ="bray", # Distance method to use, required by NMDS
             sample.label="SampleType",  # Labels will the the Sample type derived from the metadata 
             sample.order="SampleType"   # Samples will be reordered by the  "SampleType" factor
)


## Alpha diversity

# Creating graphs directly
set.seed(191931) # We specify the random seed to make the sampling reproducible
rare_global <-rarefy_even_depth(GlobalPatterns, rngseed=TRUE)
plot_richness (rare_global,color="SampleType" )

# Getting the alpha diversity indices
alfa_div = estimate_richness(rare_global, measures = c("Observed", "Shannon", "Simpson"))
alfa_div

#          Observed  Shannon   Simpson
# CL3          3488 6.529231 0.9946611
# CC1          3798 6.734570 0.9952866
# SV1          2814 6.463124 0.9962995
# ...           ...      ...       ...


# Extracting the metadata from the phyloseq object
my_factors = data.frame(sample_data(rare_global))
summary(my_factors)


# Testing for difference in richness according to sample type using
# a Kruskal Wallis test (Non-parametric ANOVA)
kruskal.test(alfa_div$Observed ~ my_factors$SampleType)

# Kruskal-Wallis rank sum test
#
# data:  alfa_div$Observed by my_factors$SampleType
# Kruskal-Wallis chi-squared = 21.57, df = 8, p-value = 0.005778


# Plotting richness by sample type
boxplot(alfa_div$Observed ~ my_factors$SampleType)


## Beta diversity

# Creating NMDS ordinations using Bray-Curtis distances
set.seed(191931)
ordination1 <-ordinate(rare_soil, method="NMDS", distance="bray")
nmds1 <-plot_ordination(rare_soil, ordination1, color="Treatment")
nmds1 


set.seed(191931)
rare_global = rarefy_even_depth(GlobalPatterns, rngseed=TRUE)
ordination2 = ordinate(rare_global, method="NMDS", distance="bray")
nmds2 = plot_ordination(rare_global, ordination2, color="SampleType")
nmds2


## Hypothesis testing

### PERMANOVA
## Creating a distance matrix from the OTU table 
dist_otus = phyloseq::distance(rare_global, method = "bray")

## Exporting metadata which includes the experimental design
table_metadata = data.frame(sample_data(rare_global))


## Running PERMANOVA with the exported data
set.seed(1)  
adonis(dist_otus ~ SampleType, data = table_metadata)

# Call:
#   adonis(formula = dist_otus ~ SampleType, data = table_metadata) 
#
# Permutation: free
# Number of permutations: 999
#
# Terms added sequentially (first to last)
#
#           Df SumsOfSqs MeanSqs  F.Model      R2 Pr(>F)    
# SampleType  8    8.2336 1.02920  5.6691 0.72736  0.001 ***
# Residuals  17    3.0863 0.18155         0.27264           
# Total      25   11.3199                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 


### Dispersion analysis

#Calculate the dispersion of the samples according to the factor SampleType
disp_SampleType = betadisper(dist_otus , table_metadata$SampleType)

# Now we can test the significance with using permutations 
set.seed(27)
permutest(disp_SampleType, pairwise=TRUE, permutations=1000)

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 1000
#
# Response: Distances
#           Df  Sum Sq  Mean Sq      F N.Perm   Pr(>F)   
# Groups     8 0.49019 0.061273 5.8591   1000 0.002997 **
# Residuals 17 0.17778 0.010458                          
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



## Session information <a name="p9"></a>

sessionInfo()

### Saving session info


data="a"
data2="B"

save.image("all_session.RData")
save(data, file="single_object.Rdata")


load("all_session.RData")
load("single_object.Rdata")
