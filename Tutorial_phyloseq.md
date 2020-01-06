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
7. [Unconstrained ordinations](#p7)
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

Cuando accedemos a la tabla de taxonomía podemos ver que los nombres de los niveles no son informativos así que podemos cambiarlos a algo mas útil
 
````
colnames(tax_table(mothur_data)) # Reporta los nombres de los taxones
#[1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6"

colnames(tax_table(mothur_data)) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family", "Genus")

colnames(tax_table(mothur_data))
#[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus" 
````

## Preparación de datos <a name="p4"></a>

Para la siguiente parte del tutorial vamos a trabajar con dos sets de datos mas interesante  
Podemos cargar el set que viene con *Phyloseq*

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

Antes de empezar con la filtración vamos a ver la distribución de las lecturas por muestra con un histograma. Para esto creamos una tabla con una columna llamada "sum" que muestra la cantidad de lecturas por muestra.

````
sample_sum_df <- data.frame(sum = sample_sums(soilrep)) 

# Media, máxima y mínima
min(sample_sum_df$sum)
median(sample_sum_df$sum)
max(sample_sum_df$sum)


# Histograma de lecturas por muestra
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color="Black", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts")  + ylab("Frequency")
````

 
![Imagen de objeto Phyloseq](https://carden24.github.com/images/distribucion_lecturas.png) 


### Selección y Filtración de muestras <a name="p4.1"></a>

Si es que hay muestras con pocas lecturas o controles podemos removerlas del set de datos usando el comando *prune_samples()*. Este comando puede utilizar una lista de muestras o también una lista de Verdaderos/FALSO

````
# Si quisiéramos remover del set de mothur las muestras que tengan menos de 3000 lecturas
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

Para seleccionar muestras usamos el comando *subset_sample()*

````
# Seleccionando muestras según los metadatos
global_soil <- GlobalPatterns %>%
  subset_samples(SampleType == "Soil") %>%
  prune_taxa(taxa_sums(.) > 0, .)
````


### Selección y filtración de OTUs <a name="p4.2"></a>

Para filtrar OTUs usamos el comando *prune_taxa()*

````
# Este comando remuevo OTUs que contribuyan con menos de 5 lecturas a todo el set
mothur_data3 = prune_taxa(names(which(taxa_sums(mothur_data) >= 5)), mothur_data2)
mothur_data3
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 247 taxa and 17 samples ]
#sample_data() Sample Data:       [ 17 samples by 1 sample variables ]
#tax_table()   Taxonomy Table:    [ 247 taxa by 6 taxonomic ranks ]
````

Si es que queremos separa grupo de OTUs de una taxonomía determinada usamos el comando *subset_taxa()*. Esto es útil si por ejemplo queremos visualizar por separado la distribución de un grupo taxonómico.

````
# Uso:
# objecto_nuevo = subset_taxa(objeto, Condición)

# Este comando separa todos las Actinobacteria seleccionando los OTUs cuyo Phylum es Actinobacteria
Actinobacteria = subset_taxa(GlobalPatterns, Phylum == "Actinobacteria")

# Este comando separa todos las Pseudomonas seleccionando los OTUs cuyo género es Actinobacteria
Pseudomonas = subset_taxa(GlobalPatterns, Genus == "Pseudomonas")
````

### Normalización <a name="p4.3"></a>
Es muy común convertir los datos a abundancias relativas, o usar otras transformaciones para cumplir asunciones estadísticas. Para esto utilizamos el comando transform_sample_counts(). Este comando a pesar que usa el objeto experimento solo modifica los datos de la tabla de OTU y no la tabla de metadatos. 

````
# Esta linea usa un objeto phyloseq y transforma las cuentas de lecturas en abundancia relativa de lecturas. 
rel_abun_soil =  transform_sample_counts(soilrep, function(x) {x/sum(x)}) 

# Esta linea convierte las cuentas con la función raíz cuadrada
rel_abun_sqrt_soil =  transform_sample_counts(rel_abun_soil, function(x) sqrt(x))
````

Otra alternativa es rareficar las muestras usando el comando *rarefey_even_depth()*, el cual tomara secuencias al azar para que todas las muestras tengan la misma cantidad de secuencias. 

````
rare_soil <-rarefy_even_depth(soilrep, rngseed=TRUE)
````

## Visualización de comunidades <a name="p5"></a>

Para gráficas de barra de abundancia relativa podemos utilizar el comando *plot_bar()* con objetos phyloseq. Esta función también acepta opciones para cambiar el color de relleno de las barras(e.g. fill="Phylum"), y puede separar las barras según otros factores 
Con este comando podemos graficar directamente o crear un objeto de tipo ggplot2 que podemos modificar.

````
# Gráfica a nivel de phyla
# Aglomeramos los datos al nivel taxonómico de Phylum
global_phyla = tax_glom(physeq = GlobalPatterns, taxrank = "Phylum")
# Convertimos a abundancia relativa
global_phyla_rel_abu = transform_sample_counts(physeq = global_phyla, function(x){x/sum(x)})

plot_bar(global_phyla_rel_abu, fill="Phylum")
````



Podemos filtrar grupos raros para que la gráfica sea mas simple

````
# Esta linea remueve grupos con abundancia menor al 5% 
global_phyla_rel_abu_clean = prune_taxa(taxa_sums(global_phyla_rel_abu) > 0.05, global_phyla_rel_abu)

phyla_plot = plot_bar(global_phyla_rel_abu_clean, fill="Phylum")
phyla_plot 
phyla_plot + facet_wrap(~SampleType, nrow=1, scales="free_x")
````

Para visualizar otros niveles de taxonomía tenemos que seleccionar subsets primero y luego graficarlos con *plot_bar()*.

````
Acidos = subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
Acidos2 = tax_glom(Acidos, taxrank = "Class")
Acidos_rel_abu = transform_sample_counts(Acidos2, function(x){x/sum(x)})
Acido_plot = plot_bar(Acidos_rel_abu, fill="Class")
Acido_plot
````
Podemos visualizar los mismos datos en un heatmap. *Phyloseq* usa ordenaciones para organizar las muestras y OTUs en vez de usar agrupamiento jerárquico como otros métodos tradicionales de heatmap.
````
# Sin reordenación de columnas
plot_heatmap(Acidos_rel_abu,  method=NULL, sample.label="SampleType", sample.order="SampleType")

# Ordenando por NMDS con distancias de Bray-Curtis
plot_heatmap(Acidos_rel_abu,  # Objeto a graficar
             method="NMDS",   # Método para crear ordenación
             distance ="bray", # Distancia para crear ordenación ya que NMDS la requiere
             sample.label="SampleType",  # Las muestras utilizan los nombre que les corresponden en la tabla de metadatos
             sample.order="SampleType"   # Las muestras están re-ordenadas según el factor "SampleType" en la tabla de metadatos
             )
````


## Diversidad alfa <a name="p6"></a>
Podemos hacer análisis de diversidad alfa directamente con la función *plot_richness()* y calcular diferentes indices de diversidad alfa con *estimate_richness()*.

````
# Graficando directamente
set.seed(191931)
rare_global <-rarefy_even_depth(GlobalPatterns, rngseed=TRUE)
plot_richness (rare_global,color="SampleType" )
plot_richness (mothur_data, color="dpw") 

# Para hacer pruebas estadísticas necesitamos extraer los estimados de diversidad alfa y los metadatos
alfa_div = estimate_richness(rare_global, measures = c("Observed", "Shannon", "Simpson"))
alfa_div

# Extraemos los metadatos
factores = data.frame(sample_data(rare_global))
summary(factores)

# Hacemos una prueba de Kruskall Wallis (ANOVA no paramétrico)
# Probamos diferencia de diversidad observada (riqueza) por tratamiento
kruskal.test(alfa_div$Observed ~ factores$SampleType)

# Graficamos la riqueza por tipo de muestras
boxplot(alfa_div$Observed ~ factores$SampleType)
````

## Ordenaciones libres <a name="p7"></a>

Para hacer ordenaciones usamos el comando *ordinate()* y luego graficamos la ordenación con *plot_ordination()*. El primer comando usar por defecto el método DCA pero puede usar también CCA, RDA, CAP, DPCoA, NMDS, MDS, y PCoA. 


````
ordination1 <-ordinate(rare_soil, method="NMDS", distance="bray")
nmds1 <-plot_ordination(rare_soil, ordination1, color="Treatment")
nmds1 

set.seed(191931)
rare_global = rarefy_even_depth(GlobalPatterns, rngseed=TRUE)
ordination2 = ordinate(rare_global, method="NMDS", distance="bray")
nmds2 = plot_ordination(rare_global, ordination2, color="SampleType")
nmds2

# Este requiere la versión mas actual de R y Phyloseq, puede que no funcione (a mi no me funciono :( )
ordination3 <-ordinate(rare_global, method="PCoA", distance="unifrac", weighted=TRUE)
w_unifrac = plot_ordination(rare_global, ordination3, color="SampleType")
w_unifrac

````

## Pruebas de hipótesis <a name="p8"></a>

### PERMANOVA <a name="p8.1"></a>
Para el análisis de variación multivariado vamos a utilizar el paquete *vegan* por lo que tenemos que exportar la tabla de OTUs y la tabla de metadatos

````
# PERMANOVA usa permutaciones al azar,
# así que si establecemos  el valor de la semilla inicial
# para que el resultado sea reproducible.
set.seed(1)  
  
# Exportación de la tabla de disimilitud de OTUs
dist_otus = phyloseq::distance(rare_global, method = "bray")

# Exportacion de la tabla de metadatos o del dise~o experimental
tabla_metadatos = data.frame(sample_data(rare_global))

summary(tabla_metadatos)

adonis(dist_otus ~ SampleType, data = tabla_metadatos)
#
#Call:
#adonis(formula = tabla_otus ~ SampleType, data = tabla_metadatos) 
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
La interpretación es que los centroides para los diferente tipo de suelos son significativamente diferente (p<0.001). 
El factor *SampleType* explica es 67.5% de la de los perfiles de OTUs.  
Si hay mas de un factor hay que incluirlos en la formula de manera secuencial:
>adonis(tabla_otus ~ Factor1 + Factor2, data = tabla_metadatos)  

### Análisis de dispersión <a name="p8.2"></a>

Luego de hacer el PERMANOVA es necesario hace un análisis de dispersión para saber si la diferencia de centroides que identifico PERMANOVA es influenciada por diferencias en dispersión.

````
# Calculamos la dispersión de la muestras según su tipo
disp_SampleType = betadisper(dist_otus , tabla_metadatos$SampleType)

# Ahora hacemos un la prueba de significancia (similar a un ANOVA) usando permutaciones para evaluar si la varianza difiere entre grupos
permutest(disp_SampleType, pairwise=TRUE, permutations=1000)
````

Los resultados sugieren que si hay diferencias en dispersión lo cual hace los resultados mas difíciles de interpretar. Se espera que esto ocurra con este set de datos porque este incluye muchos tipos de muestras incluyendo suelos (que son muy diversos) y muestras intestinales (que no los son).

### ANOSIM <a name="p8.3"></a>

El análisis de similaridad (ANOSIM) es útil cuando los grupos a comparar tienen tamaños muy diferentes como alternativa a PERMANOVA.
Este test no asume igualdad de varianza, pero es limitado en el numero de variables a evaluar y solo funciona con variables categóricas.
````
anosim_otus = anosim(dist_otus, tabla_metadatos$SampleType, permutations = 1000)
anosim_otus
plot(anosim_otus)
````

## Ordenaciones restringidas <a name="p9"></a>

Para crear ordenaciones restringidas podemos usar el comando *ordinate()*.  Vamos a trabajar con el set llamado *soilrep*. Este set de datos proviene de suelos de pastura que fueron cortados o sin cortar, con y sin calentamiento. Tiene ademas un alto nivel de replicación (24 muestras * 4 tratamientos). 


````
data(soilrep)

set.seed(191931) # Establecemos el valor inicial para el submuestreo al azar
soil2 = rarefy_even_depth(soilrep, sample.size = 1000) # Seleccionamos 1000 secuencias por muestra

# Extraemos los componentes
tabla_otus2 = t(otu_table(soil2))
dist_otus2 = phyloseq::distance(soil2, method = "bray")
tabla_metadatos2 = (data.frame(sample_data(soil2)))

# Creamos modelos estadisticos, uno con solo los tratamiento, el segundo considerando las replicas
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

# Probamos el nivel de significancia de la ordenación
anova(cap_ord1)
anova(cap_ord2)
````

Para graficar los modelos podemos usar *plot_ordination()* de nuevo.

````
# Graficamos los modelos
cap_plot1 <- plot_ordination(
  physeq = soil2, 
  ordination = cap_ord1, 
  color = "clipped",
  shape="warmed",
  axes = c(1,2)
)
cap_plot1
````

Para un gráfico mas detallada podemos obtener los datos de los centroides y agregarlos manualmente.

````
# Agregamos las variables ambientales como flechas
arrowmat <- vegan::scores(cap_ord1, display = "bp")

# Agregamos las etiquetas y creamos una tabla
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Definimos la coordenadas para agregar las flechas
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

# Definimos las etiquetas para agregar a las flechas
label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Creamos una nueva gráfica
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

## Datos de la sesión <a name="p10"></a>


Finalmente se recomienda grabar los datos de la sesión, esto es útil para poder recrear los datos completamente, especialmente si es que posteriormente se usan otras versiones de R o de alguno de los paquetes. Esta información se puede guardar al final del script.

```
sessionInfo()
```
 
Alternativamente se puede grabar todos objetos de la sesión usando 

````
save.image("Toda_la_session.RData") # Graba todos los objetos de la sesión a un archivo de extensión *.Rdata*
save(data, file="un_objeto.Rdata")  # Graba el objeto data a un archivo de extensión *.Rdata*

load("Toda_la_session.RData") # Carga todos los objetos contenidos en el archivo de extensión *.Rdata*
load("un_objeto.Rdata") # Carga todos los objetos contenidos en el archivo de extensión *.Rdata*

````

````
#Enterotipos

data(enterotype)
set.seed(20180319)

enterotype2 = prune_samples(!is.na(sample_data(enterotype)$Enterotype), enterotype)
enterotype2

dist1 = phyloseq::distance(enterotype2, method = "bray")
metadatos = sample_data(enterotype2)
summary(metadatos)

#capscale(dist1 ~ Enterotype, data=metadatos, sqrt=T)

cap_ord2 <- ordinate(
  physeq = enterotype2, 
  method = "CAP",
  distance = dist1,
  formula = ~ SeqTech + Enterotype 
)

anova(cap_ord2)
````



