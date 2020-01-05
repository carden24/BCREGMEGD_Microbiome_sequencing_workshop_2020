# Repository for
## BCREGMED Microbiome sequencing workshop 2020

## By *Erick Cardenas Poire*
-----------
This repository contains:

1. Slides
- Introduction to diversity analysis using 16S rRNA gene sequencing

2. Tutorial  for analysis of Miseq data using *Mothur*:
- We will use a modified version of a tutorial created by *Mothur* developers [(accesible here).](<https://github.com/carden24/BCREGMEGD_Microbiome_sequencing_workshop_2020/blob/master/Tutorial_mothur_Miseq_SOP.md>)
The main difference is that we will not do a phylotype analysis and use fewer samples to make the tutorial faster.  
The standard operational protocol (SOP) for *Mothur* using Miseq data can be found [here.](<http://www.mothur.org/wiki/MiSeq_SOP>)
- We require to have *Mothur* installed. The installer can be found [here](<http://www.mothur.org/wiki/Download_mothur>) (be sure to select the one for your operational system.
- The reference alignment for the 16S rRNA gene created by *Silva* can be found [here.](<http://www.mothur.org/w/images/9/98/Silva.bacteria.zip>)
- The version .9 of the training set created by the RDP for *Mothur* can be found [here.](<http://www.mothur.org/w/images/5/59/Trainset9_032012.pds.zip>)
- The files required can be found [here.](https://github.com/carden24/BCREGMEGD_Microbiome_sequencing_workshop_2020/tree/master/MiSeq_SOP_files)
- For 454 type data, *Mothur* also has a  recommended protocol [here.](<http://www.mothur.org/wiki/454_SOP>)

3. Tutorial for manipulation of *Mothur* files using *Phyloseq*
- We will use R and the *Phyloseq* package. The tutorial can be found [here](<https://github.com/carden24/BCREGMEGD_Microbiome_sequencing_workshop_2020/blob/master/Tutorial_phyloseq.md>). The microbiome data comes from the previous tutorial.
- The tutorial is based partialy on a *Phyloseq*  tutorial created by Michelle Berry found [here](http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html).

4. Other files (in the folder *Resources*):  
Publication, books, and other tools organized by themes. Including
- *Tree diversity analysis*, by Roeland Kindt and Richard Coe.  
  A short book on diversity analysis. Contains detailed explanation and code to make diversity analysis in R.
  This book can also be found freely [here.](<http://www.worldagroforestry.org/output/tree-diversity-analysis>)
- *The madness of microbiome : Attempting to find consensus “best practice” for 16S microbiome studies*, by Jolinda Pollock, Laura Glendinning, Trong Wisedchanwe, and Mick Watson.  
  A review on topics to considere when doing microbiome studies. Includes topics such as data storgae, DNA extraction, analysis, etc.  
- *Multivariate analyses in microbial ecology* by Alban Ramette. A minireview also found [here.](<http://dx.doi.org/10.1111/j.1574-6941.2007.00375.x>)
