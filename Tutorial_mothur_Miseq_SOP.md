# Microbiome data processing using *Mothur*
### By Erick Cárdenas Poiré


---

# Table of contents
1. [Introduction](#p1)
2. [Preparation](#p2)
    1. [*Mothur* installation](#p2.1)
    2. [*Mothur* reference files](#p2.2)
    3. [Sequences and related files](#p2.3)
3. [Initial processing](#p3)
   1. [Joining pairs](#p3.1)
   2. [Removal of abnormal sequences](#p3.2)
   3. [Dereplication](#p3.3)
4. [Alignment](#p4)
5. [Error removal](#p5)
   1. [Noise removal](#p5.1)
   2. [Chimera removal](#p5.2)
   3. [Unwanted lineages removal](#p5.3)
6. [OTU table creation(#p6)
7. [Phylogenetic trees creation](#p7)
8. [Data export](#p8)

## Introduction  <a name="p1"></a>

This tutorial uses the *Mothur* program to process 16S rRNA gene sequencing data. This tutorial is based on the standard operational protocol developed by Patrick Schloss which is available [here.](http://www.mothur.org/wiki/MiSeq_SOP)


**Requirements:**
- Mothur installed (see below)
- Mothur reference files (found in the MiSeq_SOP_files folder)
- Sequences and related files (found in the MiSeq_SOP_files folder)

## Preparation <a name="p2"></a>

The first step is to decide where to store all the files for the tutorial, a working directory. We want to use folder easily accesible, e.g. a folder called *mothur* in the hard drive.


### *Mothur* installation <a name="p2.1"></a>

The most recent version of *Mothur* can be found [here](https://github.com/mothur/mothur/releases/latest).  
To run this tutorial we need to download the file and uncompress it inside the working directory.

This version of *Mothur* does not have a graphic interfase. Instead uses commmands that are typed at the *mothur* terminal.
After uncompressing the file you should see something like this 

![Imagen de Mothur prompt](https://carden24.github.com/images/Mothur2.jpg) 

If that does not work, we can access windows terminal directly and initialized *Mothur* from ther.

To acess the Windows terminal presss  **[Windows] + R**  
Then type **cmd** and press **Enter**

We will be located by default at the folder of the current user ````c:\Users\erick````
To move to another folder we use the command **cd**

````
cd ..             Change to the upper directory
cd ..
cd c:\mothur      Change directly to the folder "mothur" in c:/

````

After we get to the folder where the programs are stored we can run the program by typing:
````
mothur.exe
````

In this new terminal we have to type the *Mothur* commands. *Mothur* keeps a list of the latest files created so it is not necessary to retype the whole file nane, we can just replace the name with "current".  
To know which files are the current ones we can use the commmand *get.current()*. To changes this parameters we use the command *set.current*.

*Mothur* can also use the commands directly from the Windows terminal or run a series of commands writen in a text file (one command per line).


````
# Do not run these examples
## To run one command
mothur "#fastq.info(fastq=test.fastq);get.current()"

## To run all the commands in a file
mothur stability.batch
````
The *stability.batch* has all the commands needed for the standard protocols and can be reused to analyze other files as longs as we change the content of *stability.files* which has a list of the sequencing files.


```
For linux:
mkdir mothur
cd mothur
mothur
```

### *Mothur* reference files <a name="p2.2"></a>


*Mothur* need some extra files to properly work.  These files can be found in the folder MiSeq_SOP_files.
We need to get these files and uncompress them in the working folder.

These are the files we need:

- Reference alginment from *Silva* version 102 (Silva.bacteria.zip). More recent version of the aligment may be found [here.](https://www.mothur.org/wiki/Silva_reference_files)
- Ribosomal database project (RDP) taxonomic scheme formatted for *Mothur* (Trainset14_032015.pds.zip). *Mothur* keeps a list of other taxonomics [here](https://www.mothur.org/wiki/RDP_reference_files). *Mothur* recommends using the *Greengenes* taxonomy. *Mothur* also keeps a list of more recent versions [here.] https://www.mothur.org/wiki/Greengenes-formatted_databases)


### Sequences and related files  <a name="p2.3"></a>

Finaly, we need to get the Miseq sequences and other files related to the experiment. The full tutorial uses 21 samples (two files per samples), we will use files for only nine samples to make the tutorial faster..

- Sequences in fastq format (Miseq SOP.zip)
- Experimental design file (mouse.time.design)
- Metadata (mouse.dpw.metadata)
- stability.files (a list of the sample files)

## Initial processing <a name="p3"></a>

### Joining pairs <a name="p3.1"></a>

The command *make.contigs* joins the two files that correspond to each sample as they correspond to the same DNA molecule.
The primers used to amplify the V4 region, 515F and 806R, create a product of roughly 290 bases. Since we sequence from each side with 250 bases there is a big overlap between the reads. *Mothur* joins tge reads and uses the quality information to asign bases in the overlapping regions.  


```
make.contigs(file=stability.files, processors=2)
```

The commenda created some files that are later used:
  
- *stability.trim.contigs.fasta* : The newly joined sequence
- *stability.contigs.groups* : The sample or group that each sequence belongs
- *stability.contigs.report* : A report for the concatenation process
  
To create stats for the sequences at this stage we used the command *summary.seqs*


```
summary.seqs()
 

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       249     249     0       3       1
2.5%-tile:      1       252     252     0       3       1550
25%-tile:       1       252     252     0       4       15492
Median:         1       252     252     0       4       30983
75%-tile:       1       253     253     0       5       46474
97.5%-tile:     1       253     253     6       6       60416
Maximum:        1       502     502     249     243     61965
Mean:   1       252     252     0       4
# of Seqs:      61965
```

The result shows that the majority of the sequences have between 249 and 253 bases, and that there are few sequences of up to 502 bases which suggest thtat they were not properly joined since we expected a 290 bases product.

 
### Removal of abnormal sequences <a name="p3.2"></a>

We will use the *screen.seqs* command to remove sequences according the their size and the number of ambigous bases per sequence (Ns). This step removes erronous sequences and artifacts. In this case we will sequences with any ambigous base and sequences longer than 275 bp. These parameters depend on the size of the amplified region (defined by the primers).

 
```
screen.seqs(fasta=current, group=current, summary=current, maxambig=0, maxlength=275)


>summary.seqs()
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       250     250     0       3       1
2.5%-tile:      1       252     252     0       3       1313
25%-tile:       1       252     252     0       4       13127
Median:         1       252     252     0       4       26253
75%-tile:       1       253     253     0       5       39379
97.5%-tile:     1       253     253     0       6       51192
Maximum:        1       259     259     0       11      52504
Mean:   1       252     252     0       4
# of Seqs:      52504
```
 

### Dereplication <a name="p3.3"></a>

In this step we remove identical sequences to reduce the computational load. *Mothur* keeps track of the sequences removed when creating later results. This step does not change the final results or the quality of the sequences but reduced the processing time and the memory requirements.

```
unique.seqs(fasta=current)
``` 

The protocol requires that we use the *count.seqs* command to create a table that records the identical sequences.

```
count.seqs(name=current, group=current)
```

## Alignment <a name="p4"></a> 

In this step, we aling our sequences against the reference alignment from Silva, a high-quality manualy curated alignment. This is a large files with almost 15000 sequences and over 50000 positions. To make our protocol easier, we can edit the alignment by trimming it to our region of interest (V4) and we later use this trimmed alignment with our sequences. After alignment, we will edit the new aligment which also includes our sequences.  

To edit the alignment we can use *pcr.seqs* to trim the alignment to our region of interest, then we will rename the result with *rename.file*. If we do not know the coordinates of the region for the master alignment, we can just skip this part and go directly to the alignment step with the full master alignment.

```
Do not run these commmands, we already have the final file. 
pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)
rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)

```
We can now align the sequences to our V4-trimmed master aligment from Silva. 


```
align.seqs(fasta=current, reference=silva.v4.fasta)
```

After alignment, we can run again *summary.seqs()* to create statistics for the aligned sequences.

````
summary.seqs(count=current)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1250    11546   250     0       3       1
2.5%-tile:      1968    11550   252     0       3       1313
25%-tile:       1968    11550   252     0       4       13127
Median:         1968    11550   252     0       4       26253
75%-tile:       1968    11550   253     0       5       39379
97.5%-tile:     1968    11550   253     0       6       51192
Maximum:        1982    11553   259     0       11      52504
Mean:   1967    11549   252     0       4
# of Seqs:      52504
````

These results show that the majority of the sequences start at position 1968 and end at position 11550. 
This is useful to detect sequences with too many insertions or those that start very late. We can also see that there are few sequences with many homopolymers which tend to be erroneous. 

We use this information for the *screen.seqs* command.

```
screen.seqs(fasta=current, count=current, start=1968, end=11550, maxhomop=8)

>summary.seqs(count=current)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1965    11550   250     0       3       1
2.5%-tile:      1968    11550   252     0       3       3217
25%-tile:       1968    11550   252     0       4       32164
Median:         1968    11550   252     0       4       64328
75%-tile:       1968    11550   253     0       5       96492
97.5%-tile:     1968    11550   253     0       6       125439
Maximum:        1968    13400   270     0       8       128655
Mean:   1968    11550   252.463 0       4.36666
# of unique seqs:       16298
total # of seqs:        128655

```

Now we can edit the alignment to remove position that contain only gaps because they have no information. The *screen.seqs* command also need to know which character is used in the aligment to indicate that there is no information (trump character). In our case, Silva uses ".".

```
filter.seqs(fasta=current, vertical=T, trump=.)
```

We run *unique.seqs* once again to since there are now identical sequences after the alignment and trimming step . 

```
unique.seqs(fasta=current, count=current)
```

Finaly, we run *summary.seqs()* and see that the alignment has fewer position which makes the later stepts easier. 

````
summary.seqs(count=current)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       376     249     0       3       1
2.5%-tile:      1       376     252     0       3       3217
25%-tile:       1       376     252     0       4       32164
Median:         1       376     252     0       4       64328
75%-tile:       1       376     253     0       5       96492
97.5%-tile:     1       376     253     0       6       125439
Maximum:        1       376     256     0       8       128655
Mean:   1       376     252.462 0       4.36666
# of unique seqs:       16295
total # of seqs:        128655
````

## Error removal <a name="p5"></a>

These steps improve the data quality by removing sequencing errors, host sequences, and few PCR artifacts.  

### Noise reduction <a name="p5.1"></a>

The first step is to group highly similar sequences with *pre.cluster*. This command groups sequences that differ only in few positions, in thise case no more than two differences (roughly 1 for each 100 bases). For each of the highly-similar group of sequences, the algorithm picks the most abundant sequences as the correct one and asumes that the rest are different due to sequencing errors. This step does not remove sequences.

```
pre.cluster(fasta=current, count=current, diffs=2)
```

### Chimera removal <a name="p5.2"></a>

Chimeras are artifacts created during PCR when sequences derived from one species get align with those from other specie. The resulting product does not represent the real diversity of the community and needs to be removed.
The first step detects the chimeras, and the second step removes them. The commands need to be used one after the other.

```
chimera.vsearch(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)

>summary.seqs(count=current)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       376     249     0       3       1
2.5%-tile:      1       376     252     0       3       2954
25%-tile:       1       376     252     0       4       29537
Median:         1       376     252     0       4       59073
75%-tile:       1       376     253     0       5       88609
97.5%-tile:     1       376     253     0       6       115191
Maximum:        1       376     256     0       8       118144
Mean:   1       376     252.464 0       4.37541
# of unique seqs:       2279
total # of seqs:        118144

```
vsearch removed 8.2% of the sequences (n=10511) and also removed the number of unique sequences which makes the processing easier.


### Unwanted lineages removal <a name="p5.3"></a>

In few case, the PCR can amplify sequences from the host organelles (mitochondria and chloroplast), sequences from archaeas or eukaryota(with smaller specificity and sensitivity) and other non-specific products. These sequences need to be identified and removed because they do not represent the microbial community.  

The first step is to classify all the sequences with *classify.sequences* and later remove unwanted lineages with *remove.lineage*. The first step uses a Bayesian classification method and the RDP reference taxonomy. This taxonomy is useful for the tutorial but we recommend using the *Greengenes* taxonomy which classifies up to the species level. The files required for using the *Greengenes* taxonomy with *Mothur* can be found [here.](https://www.mothur.org/wiki/Greengenes-formatted_databases)


```
classify.seqs(fasta=current, count=current, reference=trainset14_032015.pds.fasta, taxonomy=trainset14_032015.pds.tax, cutoff=80)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-mitochondria-unknown-Archaea-Eukaryota)


If working with Greengenes:
classify.seqs(fasta=current, count=current, reference=gg_13_8_99.fasta, taxonomy=gg_13_8_99.gg.tax, cutoff=80)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-mitochondria-unknown-Archaea-Eukaryota)


>summary.seqs(count=current)
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       376     249     0       3       1
2.5%-tile:      1       376     252     0       3       2950
25%-tile:       1       376     252     0       4       29496
Median:         1       376     252     0       4       58992
75%-tile:       1       376     253     0       5       88487
97.5%-tile:     1       376     253     0       6       115033
Maximum:        1       376     256     0       6       117982
Mean:   1       376     252.465 0       4.37191
# of unique seqs:       2259
total # of seqs:        117982

```
The original Miseq tutorial explains how to use mock communities (artificial mixtures with known abundances) to calculate error rates. If you plan to use these type of positive controls  you can read the protocols [tutorial original de Mothur](https://www.mothur.org/wiki/MiSeq_SOP). 
In this case we will simply remove the mock sequences from the analysis using the *remove.seqs* command.


```
remove.groups(count=current, fasta=current, taxonomy=current, groups=Mock)
```

## OTU table creation <a name="p6"></a>

The next objetive is to group sequences in operational taxonomic units (OTUs), sequence groups defined by the similarity among them. Usually, the first step is to compare all sequences against each other to create a distance matrix (using *dist.seqs*) and later use the matrix to create group of sequences (using *cluster*). This method is slow when working with thousand of sequences since the matrix uses too much memory. 

````
Do not run these
dist.seqs(fasta=current, cutoff=0.03)
cluster(column=current, count=current)

````
A more appropriate altenative is to use *cluster.split*. This splits the matrix into groups defined by their taxonomic affiliation and later use the groups for clustering. We will use the taxonomic classification results to create group at taxonomic level of order  and later use the *opticlust* method to create OTUs. 

The *cluster.split* command needs to know the similarity threshold for the creation of the distance matrix. We will use "cutoff=0.03", which is appropriate for this method. For the *average neighbour* methodm, a higher parameter of 0.15 is recommended for working at species level.


```
cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.03)
```

Now we can use the distance matrices to generate the OTUs table. In this case, we want to work at species level which is approximately 3%. 

```
make.shared(list=current, count=current, label=0.03)
```
The generated table is all we need for diversity analysis. 

The last step in the protocol is to classify each OTU. We will use the *classify.otu* command to generate a consensus classification for each OTU. 

```
classify.otu(list=current, count=current, taxonomy=current, label=0.03)
```

## Phylogenetic trees creation <a name="p7"></a>

Some methods, such as UNIFRAC, require to know the place of each OTU in this project in a phylogenetic tree. To create this tree we use a distance-based tree method using our alignment and its corresponding distance matrix (Neighbour-joining method).

```
dist.seqs(fasta=current, output=lt)
clearcut(phylip=current)
```

## Data export <a name="p8"></a>

*Mothur* generates many intermediate files with complicated names. An easy way to export the key information is to create a biom file. The biom format stores the information in a standard way which is easy to process by many other programs. More information on this format can be found [here.](http://biom-format.org/documentation/biom_format.html)


````
make.biom(shared=current, constaxonomy=current, metadata=mouse.dpw.metadata)
````  


## Phylotype analysis   [Opcional]

If we are not interested in creating an OTU table (species defined by their similarity), but work with species defined by their taxonomic classification (phylotypes), we can use the following commands on processed data (without noise and artifacts).

````
phylotype(taxonomy=current, label=1);
make.shared(list=current, count=current);
classify.otu(list=current, taxonomy=current);
````

The first command group sequences according to their taxonomy, the option "label=" specifies that we will use the most specific taxonomy. If we do not specy this, the results are going to be list of sequences at different levels of taxonomy.  
The second command, creates the phylotype x samples table. The last command creates another table with the consensus classification for each phylotype. We need to be careful to select the correct files since the names are very similar.

```
"Shared" files:
# Based on OTUs
stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.shared
# Based on phylotypes
stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.shared           <-Basado en filotipos

# Consensus taxonomy files:
# Based on OTUs

stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.1.cons.tax.summary
# Based on phylotypes
stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.tax.summary

```

**Final recommendations**
Once we finish the processing we recommend saving the original fastq files (we usually need to submit them for publication), the list of commands that we used (to replicate the process), the computer logs (".logfile"), the OTU table (".shared"), the consensus classification for the OTUs (".cons.tax.summary"), the metadata, and the phylogentic tree (".phylip.tre") if we are interested in using UNIFRAC. The rest of the files can be deleted or compressed for storage.
 









