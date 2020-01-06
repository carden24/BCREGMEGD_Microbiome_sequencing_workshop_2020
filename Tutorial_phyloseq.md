# Microbial community analysis with *Phyloseq*
### By Erick Cardenas Poire

---

# Table of contents
1. [Introduction](#p1)
2. [Preparation](#p2)
3. [Creation of a Phyloseq object](#p3)
    1. [Data import](#p3.1)
    2. [Accessing components](#p3.2)
4. [Data preparation](#p4)
    1. [Sample filtering](#4.1)
    2. [OTU filtering](#p4.2)
    3. [Normalization](#p4.3)
5. [Visualization of communities](#p5)
6. [Alpha diversity](#p6)
7. [Beta diversity](#p7)
8. [Hypothesis testing](#p8)
    1. [PERMANOVA](#p8.1)
    2. [Dispersion analysis](#p8.2)
    3. [ANOSIM](#p8.3)
9. [Constrained analysis](#p9)
10. [Session information](#p10)

## Introduction <a name="p1"></a>

This tutorial covers the analysis of microbial community data using package *Phyloseq* in R. The tutorial is partialy based on the *Phyloseq* tutorial developed by Michelle Berry available [here.](http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html)


*Phyloseq* is an R package focused on analysis of microbiome census data. The original paper can be found in the *Resources* folder in this repository. For more information on *Phyloseq* go [here.](http://joey711.github.io/phyloseq/)

The main advantage of using *Phyloseq* it that it creates an object that groups and organizes all the key data for a microbiome study: the OTUs table, metadata, taxonomic classification of the OTUs, and a phylogenetic tree for the OTU representative sequences. *Phyloseq* provides specialized functions for import, export, manipulation, and visualization of microbiome data.

**Phyloseq-class object scheme**  
![Imagen de objeto Phyloseq](https://carden24.github.com/images/Phyloseq.jpg)  


## Preparation <a name="p2"></a>

The first recommended step is to set the working directory. This is the directory where *R* looks for files and where it sabes the exported results. We recommend that you create a working directory for each major analysis and that you also store there the analysis scripts.

Open *RStudio*, and create a new *R* script with **Ctrl** + **Shift** + **N**. This is where we will write the commands and execute each line with **Ctrl** + **R**. In *R* the symbol **#** is used to add comments to that line of code. Everything that comes after the "#" is ignored by the program but it is useful for us to know what we did for each line or block code.


````
setwd("E:/Libraries/Dropbox/tutorial/") # This is my working directory, change it to something more useful for you
getwd() # This command report me with current working directory

````

The second step is to load the *Phyloseq* library to access the functions required for the tutorial.

````
library(phyloseq)
library(ggplot2) # This package is required for data visualization
library(vegan) # This package is required for the diversity analysis
library(tidyr) # These two package help with data manipulation
library(dplyr))

````

Finaly, we will load some custom functions that are stored in the "miseqR.R". These functions come from the original tutorial and are just useful tools. The file should be saved in the working directory for *R* to find it. The files required for this tutorial can be found in the "Phyloseq_files" folder of this[workshop repository.](https://github.com/carden24/BCREGMEGD_Microbiome_sequencing_workshop_2020)


````
source("miseqR.R")
````

## Creation of a Phyloseq object <a name="p3"></a>

### Data import <a name="p3.1"></a>

*Phyloseq* has a handy function to import files from *Mothur*. So, we need to move the final *Mothur* files to the *R* working directoy and use the code below. The files can also be found in the "Phyloseq_files" folder [here.] (https://github.com/carden24/BCREGMEGD_Microbiome_sequencing_workshop_2020/tree/master/Phyloseq_files)


````
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
# otu_table()   OTU Table:         [ 522 taxa and 19 samples ]
# sample_data() Sample Data:       [ 19 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 522 taxa by 6 taxonomic ranks ]
````

Alternatively we can add each element individually if we do not have the biom file. Please note that the biom file did not have the phylogenetic tree.

````
# Importing Mothur files
mothur_data0 <- import_mothur(mothur_shared_file = shared_file,
                             mothur_constaxonomy_file = tax_file)

# Adding metadata

metadata = read.delim(metadata_file, header=T, row.names=1)
metadata = sample_data(metadata) # Converts the metadata table into a metadata table ready for *Phyloseq*
mothur_data = merge_phyloseq(mothur_data0, metadata)

mothur_data
# otu_table()   OTU Table:         [ 522 taxa and 19 samples ]
# sample_data() Sample Data:       [ 19 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 522 taxa by 6 taxonomic ranks ]
````

### Accessing components <a name="p3.2"></a>

The *mothur_data* object is a *Phyloseq*-class object of experiment-level. To see the structure of this object we simply type its name and run the line. 


````
mothur_data
# otu_table()   OTU Table:         [ 522 taxa and 19 samples ]
# sample_data() Sample Data:       [ 19 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 522 taxa by 6 taxonomic ranks ]
```

Once we have this object we can use several function to access the different components of the object.

```
# Basic functions, just put the name of the object inside the parenthesis
 
otu_table(mothur_data)			# Reports the OTU table
sample_data(mothur_data) 		# Reports the sample information such as experimental design or metadata
tax_table(mothur_data)   		# Reports the taxonomic classification of each OTUs

ntaxa(mothur_data) 			# Number of OTUs
nsamples(mothur_data) 			# Number of samples

sample_names(mothur_data)			# Sample names
taxa_names(mothur_data) 			# Taxa names

taxa_sums(mothur_data)  	# Sum of abundances or an OTU for all the samples
````

When we access the taxonomy table we can see that the names of the ranks are not information so we are going to change them to something more useful.
 
````
colnames(tax_table(mothur_data)) # Reports the names of the taxonomic ranks
#[1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6"

# Replacing the names
colnames(tax_table(mothur_data)) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family", "Genus")

colnames(tax_table(mothur_data))
#[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus" 
````

## Data preparation <a name="p4"></a>

Podemos cargar el set que viene con *Phyloseq*
For the next part of the tutorial we will use more interesting data that come as example in the *Phyloseq* package. We will load the "GlobalPatterns" and "soilrep" datasets. The first one comes from a 2011 [publication] (http://www.pnas.org/content/108/suppl.1/4516.short) comparing 25 environmental samples and three known “mock communities” – a total of 9 sample types. The second dataset come from a 2011 [paper](http://www.nature.com/ismej/journal/v5/n8/full/ismej201111a.html) which compared 24 separate soil microbial communities under four treatment conditions.

````
data("GlobalPatterns")
GlobalPatterns
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
#sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]

data("soilrep")
soilrep
#otu_table()   OTU Table:         [ 16825 taxa and 56 samples ]
#sample_data() Sample Data:       [ 56 samples by 4 sample variables ]
````

The first thing we can do is to check the read distribution for the samples. We will create a table with a column called "sum" and visualize this data.

```
sample_sum_df <- data.frame(sum = sample_sums(soilrep)) 

# Exploring mean, maximum, and minimum
min(sample_sum_df$sum)
median(sample_sum_df$sum)
max(sample_sum_df$sum)


# Reads per samples histogram
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color="Black", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts")  + ylab("Frequency")
```

![Read distribution graph](https://carden24.github.com/images/distribucion_lecturas.png) 


### Selecting and filtering samples <a name="p4.1"></a>

If there are samples with very few sequences or controls, we can remove them from the dataset using the *prune_samples()* command. This command accepts as input a list of samples or a list with TRUE/FALSE.

````
# To remove samples with fewer than 300 reads
mothur_data2 = prune_samples(names(which(sample_sums(mothur_data) >= 3000)), mothur_data)

mothur_data
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 522 taxa and 19 samples ]
#sample_data() Sample Data:       [ 19 samples by 1 sample variables ]
#tax_table()   Taxonomy Table:    [ 522 taxa by 6 taxonomic ranks ]

mothur_data2
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 522 taxa and 17 samples ]
#sample_data() Sample Data:       [ 17 samples by 1 sample variables ]
#tax_table()   Taxonomy Table:    [ 522 taxa by 6 taxonomic ranks ]
````

To select samples we use*subset_sample()*

````
# Selecting samples according to their metadata, in this case we want only soil samples
global_soil <- GlobalPatterns %>%
  subset_samples(SampleType == "Soil") %>%
  prune_taxa(taxa_sums(.) > 0, .)
````


### Selecting and filtering OTUs <a name="p4.2"></a>

To filter OTUs we use *prune_taxa()*

````
# The following line removes OTUs that have five or fewer reads in the whole dataset
mothur_data3 = prune_taxa(names(which(taxa_sums(mothur_data) >= 5)), mothur_data2)

mothur_data3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 247 taxa and 17 samples ]
#sample_data() Sample Data:       [ 17 samples by 1 sample variables ]
#tax_table()   Taxonomy Table:    [ 247 taxa by 6 taxonomic ranks ]
````

If we want to select a group of OTUs with a specific taxonomy, we use *subset_taxa()*. This is useful if we want to visualize separately the distribution of a taxonomic group e.g. to visualize a low abundance group that will get lost in a general distribution graph.


````
# Usage:
# new_object = subset_taxa(object, Condition)

# This line will select all OTUs whose Phylum is Actinobacteria
Actinobacteria = subset_taxa(GlobalPatterns, Phylum == "Actinobacteria")

# This line will get all Pseudomonas by selecting OTUs whose Genus is Pseudomonas
Pseudomonas = subset_taxa(GlobalPatterns, Genus == "Pseudomonas")
````

### Normalization <a name="p4.3"></a>

Sometimes it is useful to convert our count data into relative abundances or to use other transformations to satisfy statistical assumptions. For this we will use *transform_sample_counts()*. This commands modifies only the OTU table and not the metadata table

````
# The next line  will take a phyloseq object and transform counts into relative abundance (0 - 1 range)
rel_abun_soil =  transform_sample_counts(soilrep, function(x) {x/sum(x)}) 

# The next line will do a square-root transformation of the counts
rel_abun_sqrt_soil =  transform_sample_counts(rel_abun_soil, function(x) sqrt(x))
````

Another alternative is to rarify the counts using *rarefey_even_depth()*, this randomly take sequences from each sample so all of them have the same number of sequences. This is useful for alpha diversity calculations if there are major differences in sequence depth.

````
rare_soil <-rarefy_even_depth(soilrep, rngseed=TRUE)
````

## Community visualization <a name="p5"></a>

A common way to visualize the whole community is by using bar plot with the *plot_bar()* command. This function also accepts options to change the bar colors or to separate the plot according to experimental factors present in the phyloseq object. Plot object can be further modified with ggplot2 parameters.


````
# Creating a phyla level plot
# Collapse data at Phylum level
global_phyla = tax_glom(physeq = GlobalPatterns, taxrank = "Phylum")
# Converting to relative data
global_phyla_rel_abu = transform_sample_counts(physeq = global_phyla, function(x){x/sum(x)})

plot_bar(global_phyla_rel_abu, fill="Phylum")
````

We could additional filter rare groups to simplify the plot

````
# This line remove groups with abundance smaller than 5% 
global_phyla_rel_abu_clean = prune_taxa(taxa_sums(global_phyla_rel_abu) > 0.05, global_phyla_rel_abu)

phyla_plot = plot_bar(global_phyla_rel_abu_clean, fill="Phylum")
phyla_plot 
phyla_plot + facet_wrap(~SampleType, nrow=1, scales="free_x")
````

To visualize some parts of the communities we need to subset them first and then use *plot_bar*

````
Acidos = subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
Acidos2 = tax_glom(Acidos, taxrank = "Class")
Acidos_rel_abu = transform_sample_counts(Acidos2, function(x){x/sum(x)})
Acido_plot = plot_bar(Acidos_rel_abu, fill="Class")
Acido_plot
````

We can also visualize the same data using heatmaps via the *plot_heatmap* command. *Phyloseq* uses ordination methods to group samples instead of hierarchical clustering, commonly used in other heatmap functions.


```
# A heatmap without sample reordering
plot_heatmap(Acidos_rel_abu,  method=NULL, sample.label="SampleType", sample.order="SampleType")


# A heatmap with the samples reordered according to their sample similarity
plot_heatmap(Acidos_rel_abu,  # Object to plot
             method="NMDS",   # Method for the ordination
             distance ="bray", # Distance method to use, required by NMDS
             sample.label="SampleType",  # Labels will the the Sample type derived from the metadata 
             sample.order="SampleType"   # Samples will be reordered by the  "SampleType" factor
             )
````


## Alpha diversity <a name="p6"></a>

We can calculate and visualize alpha diversity using *estimate_richness()* (getting multiple diversity indices) and *plot_richness()*.


````
# Creating graphs directly
set.seed(191931) # We specify the random seed to make the sampling replicable
rare_global <-rarefy_even_depth(GlobalPatterns, rngseed=TRUE)
plot_richness (rare_global,color="SampleType" )
plot_richness (mothur_data, color="dpw") 

# Getting the alpha diversity indices
alfa_div = estimate_richness(rare_global, measures = c("Observed", "Shannon", "Simpson"))
alfa_div

# Extracting the metadata from the phyloseq object
my_factors = data.frame(sample_data(rare_global))
summary(my_factors)

# Testing for difference in richness according to sample type using
# a Kruskal Wallis test (Non-parametric ANOVA)
kruskal.test(alfa_div$Observed ~ my_factors$SampleType)

# Plotting richness by sample type
boxplot(alfa_div$Observed ~ my_factors$SampleType)
````

## Beta diversity <a name="p7"></a>

We use unconstrained ordinations to analyze the diversity across samples i.e. their beta diversity. 
We create ordination using *ordinate()* and then graph the ordination results with *plot_ordination()*. 
The first commands uses DCA by default but we can also use CCA, RDA, CAP, DPCoA, NMDS, MDS, and PCoA. 


````
ordination1 <-ordinate(rare_soil, method="NMDS", distance="bray")
nmds1 <-plot_ordination(rare_soil, ordination1, color="Treatment")
nmds1 

set.seed(191931)
rare_global = rarefy_even_depth(GlobalPatterns, rngseed=TRUE)
ordination2 = ordinate(rare_global, method="NMDS", distance="bray")
nmds2 = plot_ordination(rare_global, ordination2, color="SampleType")
nmds2
````

## Hypothesis testing <a name="p8"></a>

### PERMANOVA <a name="p8.1"></a>

To test the effect of a experimental factor on the similarity of samples we will do a Permutational multivariate analysis of variation (PERMANOVA) using the *vegan* package in *R*. We need to first export the OTU table and the metadata as objects.

````
# PERMANOVA uses random permutations, so we should use a specific random seed if we want to have a fully reproducible work.
set.seed(1)  
  
# Creating a distance matrix from the OTU table 
dist_otus = phyloseq::distance(rare_global, method = "bray")

# Exporting metadata which includes the experimental design
table_metadata = data.frame(sample_data(rare_global))

summary(table_metadata)

adonis(dist_otus ~ SampleType, data = table_metadata)
#
#Call:
#adonis(formula = tabla_otus ~ SampleType, data = table_metadata) 
#
#Permutation: free
#Number of permutations: 999
#
#Terms added sequentially (first to last)
#
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#SampleType  8    7.7744 0.97180  4.4259 0.67562  0.001 ***
#Residuals  17    3.7327 0.21957         0.32438           
#Total      25   11.5070                 1.00000           
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 
````
The interpretation of the results is that the centroids for the different soil types are significantly different (p<0.001). The *SampleType* factor explaines 67.5% of the OTUs profile variation.  
If there are more factor we can include them in a sequential manners:
>adonis(tabla_otus ~ Factor1 + Factor2, data = table_metadata)  

### Dispersion analysis <a name="p8.2"></a>

After a PERMANOVA is commonly used to run a dispersion analysis, to see if the differences in centroids found in PERMANOVA are influenced by difference in dispersion (how separated a group of samples are). 

````
#Calculate the dispersion of the samples according to the factor SampleType
disp_SampleType = betadisper(dist_otus , table_metadata$SampleType)

# Now we can test the significance with using permutations 
permutest(disp_SampleType, pairwise=TRUE, permutations=1000)
````

The results are significant (samples from different sample types have different dispersion), which makes the PERMANOVA results harder to interpret. We expect this in this particular set because it includes soils samples (highly variable and diverse) and gut samples (low diversity and variation).

### ANOSIM <a name="p8.3"></a>

The analysis of similarity (ANOSIM) is a useful alternative to PERMANOVA when the groups being compared have very different sizes. The test evaluates if the similarity within groups is smaller than the similarity between groups. ANOSIM does not assume equal variation but it is limited in the number of variables to evaluate and only works with categorical variables.
````
anosim_otus = anosim(dist_otus, tabla_metadatos$SampleType, permutations = 1000)
anosim_otus
plot(anosim_otus)
````

## Constrained ordinations <a name="p9"></a>

We use again *ordinate()* for constrained ordinations. Constrained ordinations try to represent the majority of the multivariate variation explained by our explanatory variables instead of all the variability present.
We will work for this examples with the *soilrep* dataset. This sets come from pasture soils that were cut/uncut and were warmed/unwarmed (four treatments in total).
. 

````
data(soilrep)

set.seed(191931) # Using a specific random seed for reproducibility
soil2 = rarefy_even_depth(soilrep, sample.size = 1000) # Randomly selecting 1000 sequences per sample

# Extracting components
table_otus2 = t(otu_table(soil2))
dist_otus2 = phyloseq::distance(soil2, method = "bray")
table_metadatos2 = (data.frame(sample_data(soil2)))

# Creating CAP models, one with treatments only, the other one also with replicates
cap_ord1 <- ordinate(
  physeq = global2, 
  method = "CAP",
  distance = dist_otus2,
  formula = ~ Treatment + warmed 
)

cap_ord2 <- ordinate(
  physeq = global2, 
  method = "CAP",
  distance = dist_otus2,
  formula = ~ Treatment + warmed + Sample 
)

# Testing significant of models
anova(cap_ord1)
anova(cap_ord2)
````

To plot the ordinations we use *plot_ordination()* again.

````
# CAP plot 1
cap_plot1 <- plot_ordination(
  physeq = soil2, 
  ordination = cap_ord1, 
  color = "clipped",
  shape="warmed",
  axes = c(1,2)
)
cap_plot1
````

We can further modify the graph by getting the centroids data from the results and adding them manually .

````
# Save environmental variables data to display them latter as arrows
arrowmat <- vegan::scores(cap_ord1, display = "bp")

# Adding labels and converting it to a data frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Creating arrow coordinates, they start at (0, 0)
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

# Creating labels to add to the arrows
label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Creating the final graph. We start with the previous plot and add text and arrows
cap_plot1 + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )
````

## Session information <a name="p10"></a>

We always recommend to save the session information (the version of R and the packages used). This is useful to fully reproduce the analysis and results since different version of R and packages can create slightly different results.  

```
sessionInfo()
```
We highly recommend saving this info (at the end of the script for example).  

We can also save a specific R object with *save.image* or all of the session objects with *save.image*


````
save.image("All_session_data.RData") # Saves all the object in the session as a *RData* file.
save(data, file="an_object.Rdata")  # Saves the "data" object from the session as a *RData* file

load("All_session_data.RData") # Loads all the object in the file with extension *.Rdata*
load("an_object.Rdata") # Same as above, but only one object is loaded since we saved only one

````

````
# Using the Enterotypes dataset

data(enterotype)
set.seed(20180319)

enterotype2 = prune_samples(!is.na(sample_data(enterotype)$Enterotype), enterotype)
enterotype2

dist1 = phyloseq::distance(enterotype2, method = "bray")
metadata_ent = sample_data(enterotype2)
summary(metadata_ent)


cap_ord2 <- ordinate(
  physeq = enterotype2, 
  method = "CAP",
  distance = dist1,
  formula = ~ SeqTech + Enterotype 
)

anova(cap_ord2)
````
