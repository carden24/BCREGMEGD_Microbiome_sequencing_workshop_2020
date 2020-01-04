  rm(list=ls()) 


setwd("E:/Libraries/Dropbox/tutorial/")



library(vegan)
library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)



install.packages("dplyr")

mothur(data)



# Crear variables para los archivos que exportamos de *Mothur*

shared_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.shared"
tax_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy"
metadata_file = "mouse.dpw.metadata"
biom_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.biom"
# Arboles importados de mothur usan los nombres de las secuencias y no de los OTUs.
# Hay que cambiar esto manualmente.
#tree_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.phylip.tre"


# Importar un objeto biom
mothur_data = import_biom(biom_file)


# Alternativamente podemos agregar cada elemento individualmente
# si que no tenemos un archivo biom

# Import mothur data
mothur_data0 <- import_mothur(mothur_shared_file = shared_file,
                             mothur_constaxonomy_file = tax_file)

# Y ahora agregamos los metadatos

metadata = read.delim(metadata_file, header=T,row.names = 1)
metadata = sample_data(metadata) # Convierta la tabla de metadatos en una tabla para phyloseq
mothur_data = merge_phyloseq(mothur_data0, metadata)

mothur_data
# otu_table()   OTU Table:         [ 522 taxa and 19 samples ]
# sample_data() Sample Data:       [ 19 samples by 1 sample variables ]
# tax_table()   Taxonomy Table:    [ 522 taxa by 6 taxonomic ranks ]


# Una vez que tenemos este objeto podemos usar diferente funciones para acceder a los diferentes componentes del bjeto

otu_table(mothur_data)   # Reporta la tabla de OTUs
sample_data(mothur_data) # Reporta la informacion sobre las muestras como metadatos o datos experimentales
tax_table(mothur_data)   # Reporta la tabla de taxonomica de los otus


ntaxa(mothur_data) # Reporta el numero de OTUs
#[1] 522

nsamples(mothur_data) # Reporta el numero de muestras
#[1] 19


sample_names(mothur_data) # Reporta los nombres de las muestras
taxa_names(mothur_data) # Reporta los nombres de los taxones


taxa_sums(mothur_data) # Reporta la suma de abundancias de un  otu para todas las muestras


# Cuando accedemos a la tabla de taxonomia podemos ver que los nombres de los niveles no son informativos asi que podemos cambiarlos a algo mas util
 

colnames(tax_table(mothur_data)) # Reporta los nombres de los taxones
# [1] "Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6"

colnames(tax_table(mothur_data)) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family", "Genus")

colnames(tax_table(mothur_data))
# [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus" 


Para la siguiente parte del tutorial vamos a trabajar con dos sets de datos mas interesante

Podemos cargar el set que viene con Phyloseq


"""
data("GlobalPatterns")
GlobalPatterns
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]


data("soilrep")
>soilrep
otu_table()   OTU Table:         [ 16825 taxa and 56 samples ]
sample_data() Sample Data:       [ 56 samples by 4 sample variables ]
"""


### Exploracion

Creamos una tabla con una columna que muestra la cantidad de lecturas por cada muestra.
y asignamos este valor a una nueva columna llamada "sum" en la tabla.


sample_sum_df <- data.frame(sum = sample_sums(soilrep)) 


# Media, maxima y minima
min(sample_sum_df$sum)
median(sample_sum_df$sum)
max(sample_sum_df$sum)


# Histograma de lecturas por muestra
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color="Black", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts")  + ylab("Frequency")


# Filtracion de muestras
mothur_data2 = prune_samples(names(which(sample_sums(mothur_data) >= 3000)), mothur_data)

mothur_data
mothur_data2

### Filtracion de OTUs

# Este comando remuevo OTUs que contribuyan con menos de 5 lecturas a todo el set
mothur_data3 = prune_taxa(names(which(taxa_sums(mothur_data) >= 5)), mothur_data2)
mothur_data3




# Abundancia a nivel de phyla

# Aglomeramos los datos al nivel taxonomico de Phylum
global_phyla = tax_glom(physeq = GlobalPatterns, taxrank = "Phylum")
# Convertimos a abundancia relativa
global_phyla_rel_abu = transform_sample_counts(physeq = global_phyla, function(x){x/sum(x)})
# Convertimos en formato largo
global_phyla_rel_abu_long = psmelt(global_phyla_rel_abu)

# Grafica a nivel de phyla
ggplot(global_phyla_rel_abu_long, aes(x = Sample, y = Abundance, fill = Phylum)) + 
    geom_bar(stat = "identity") +
  ylab("Relative Abundance") 


# Podemos filtrar grupor raros para que la grafica sea mas simple
# Esta linea remuve abundancia de grupos con menos de 5% 
global_phyla_rel_abu_long_filered = subset(global_phyla_rel_abu_long, Abundance > 0.05)

# Grafica a nivel de phyla
ggplot(global_phyla_rel_abu_long_filered, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  ylab("Relative Abundance (Phyla > 5%)") + facet_wrap(~SampleType, nrow=1, scales = "free_x")




# Una forma mas eficiente es utilizando los paquetes tidyR y dpplyr para manipular datos
global_phylum <- GlobalPatterns %>%
  tax_glom(taxrank = "Phylum") %>%                     	# Aglomera a nivel de phylum
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 	# Convierte en abudancia relativa
  psmelt() %>% 											# Convierte en formato largo`
  filter(Abundance > 0.05) %>%                         	# Filtra grupos raros , <5%
  arrange(Phylum)  # sort by Phylum name

ggplot(global_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  ylab("Relative Abundance (Phyla > 5%)") + facet_wrap(~SampleType, nrow=1, scales = "free_x")


sample_data(global_phyla)

# Seleccionando por muestroas
global_soil <- GlobalPatterns %>%
  subset_samples(SampleType == "Soil") %>%
  prune_taxa(taxa_sums(.) > 0, .)
global_soil


# Este comando separa todos las Actinobacteria 
Actinobacteria = subset_taxa(GlobalPatterns, Phylum == "Actinobacteria")
Actinobacteria



# Este comando separa todos las Pseudomonas
Pseudomonas = subset_taxa(GlobalPatterns, Genus == "Pseudomonas")
Pseudomonas



# Normalization


rel_abun_soil =  transform_sample_counts(soilrep, function(x) {x/sum(x)}) 

rel_abun_sqrt_soil =  transform_sample_counts(rel_abun_soil, function(x) sqrt(x))

plot_bar(Pseudomonas
         , fill="Species")





plot_bar(global_phyla_rel_abu, fill="Phylum")


global_phyla_rel_abu_clean = prune_taxa(taxa_sums(global_phyla_rel_abu) > 0.05, global_phyla_rel_abu)
phyla_plot = plot_bar(global_phyla_rel_abu_clean, fill="Phylum")
phyla_plot 
phyla_plot + facet_wrap(~SampleType, nrow=1, scales="free_x")



#

Acidos = subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
Acidos2 = tax_glom(Acidos, taxrank = "Class")
Acidos_rel_abu = transform_sample_counts(Acidos2, function(x){x/sum(x)})
Acido_plot = plot_bar(Acidos_rel_abu, fill="Class")
Acido_plot


#Heatmap

# Sin reordenacion de muestras
plot_heatmap(Acidos_rel_abu,  method=NULL, sample.label="SampleType")
# 
plot_heatmap(Acidos_rel_abu,  # Objeto a graficar
             method="NMDS",   # Metodo para crear ordenacion
             distance ="bray", # Distancia para crear ordenacion ya que NMDS la requiere
             sample.label="SampleType",  # Las muestras utilizan los nombre que les corresponden en la tabla de metadatos
             sample.order="SampleType"   # Las muestras estan reordenadas segun el factor "SampleTupe" en la tabla de metadatos
             )
  

# Diversidad alfa

set.seed(191931)
rare_global <-rarefy_even_depth(GlobalPatterns, rngseed=TRUE)
plot_richness (rare_global,color="SampleType" )
plot_richness (mothur_data, color="dpw") 
  
alfa_div = estimate_richness(rare_global, measures = c("Observed", "Shannon", "Simpson"))
alfa_div

# Extraemos los metadatos
factores = data.frame(sample_data(rare_global))
summary(factores)

# Hacemos una prueba de Kruskall Wallis (ANOVA no parametrico)
# Probamos diferencia de diversidad observada (riqueza) por tratamiento
kruskal.test(alfa_div$Observed ~ factores$SampleType)

boxplot(alfa_div$Observed ~ factores$SampleType)


# NMDS
  
set.seed(191931)
rare_soil <-rarefy_even_depth(soilrep, rngseed=TRUE)
ord1 <-ordinate(rare_soil, method="NMDS", distance="bray")
p1 <-plot_ordination(rare_soil, ord1, color="Treatment")
p1
  
  
set.seed(191931)
rare_global <-rarefy_even_depth(GlobalPatterns, rngseed=TRUE)
ord2 <-ordinate(rare_global, method="NMDS", distance="bray")
p2 <-plot_ordination(rare_global, ord2, color="SampleType")
p2
  
  
ordination3 <-ordinate(rare_global, method="PCoA", distance="unifrac", weighted=TRUE)
w_unifrac = plot_ordination(rare_global, ordination3, color="SampleType")
w_unifrac
  
  
  ig <- make_network(GlobalPatterns, type = "taxa", max.dist=0.5)
  plot_network(ig, GlobalPatterns, color="SampleType", line_weight=0.3, label=NULL)
  
  sample_data(GlobalPatterns)
  
  ## PERMANOVA
  
  
  # Permanova usar permutaciones al azar,
  # asi que si establecemos  el valor de la semilla inicial
  # para que el resultado sea reproducible.
  set.seed(1)  
  
  # Exportacion de la tabla de OTUs
  dist_otus = phyloseq::distance(rare_global, method = "bray")

  
  # Exportacion de la tabla de metadatos o del dise~o experimental
  tabla_metadatos = data.frame(sample_data(rare_global))

summary(tabla_metadatos)
  
adonis(tabla_otus ~ SampleType, data = tabla_metadatos)

# Calculamos la dispersion de la muestras segun su tipo
disp_SampleType = betadisper(dist_otus , tabla_metadatos$SampleType)

# Ahora hacemos un la prueba de significancia (similar a un ANOVA) usando permutaciones para evaluar si las variacias difieren por grupo
permutest(disp_SampleType, pairwise=TRUE, permutations=1000)


## ANOSIM

El analisis de similaridad (ANOSIM) es util cuando los grupos
tienen tama~os muy diferentes como alternativa a PERMANOVA
Este test no asume variancias iguales, pero es mas similato en el numero de variables a evaluar y solo funciona con variables categoricas.

anosim_otus = anosim(dist_otus, tabla_metadatos$SampleType, permutations = 1000)
anosim_otus
plot(anosim_otus)


# CAP

# Este set de datos proviene de suelos de grases que fueron cortados o sin cortar, con y sin calentamiento.
# Tiene un alto nivel de replicacion (24 muestras * 4 tratamientos). 
data(soilrep)


# Establecemos el valor inicial para el submuestreo al azar
set.seed(191931)
soil2 = rarefy_even_depth(soilrep, sample.size = 1000) # Seleccionamos 1000 sequencias por muestras


# Extraemos los componentes
tabla_otus2 = t(otu_table(soil2))
dist_otus2 = phyloseq::distance(soil2, method = "bray")
tabla_metadatos2 = (data.frame(sample_data(soil2)))

# Creamos modelos

cap_ord1 <- ordinate(
  physeq = soil2, 
  method = "CAP",
  distance = dist_otus2,
  formula = ~ Treatment + warmed 
)

cap_ord2 <- ordinate(
  physeq = soil2, 
  method = "CAP",
  distance = dist_otus2,
  formula = ~ Treatment + warmed + Sample 
)

# Probamos el nivel de significancia de la ordenacion
anova(cap_ord1)
anova(cap_ord2)



# Graficamos los modelos
cap_plot2 <- plot_ordination(
  physeq = soil2, 
  ordination = cap_ord2, 
  color = "clipped",
  shape="warmed",
  axes = c(1,2)
)
cap_plot2

# Agregamos las variables ambientales como flechas
arrowmat <- vegan::scores(cap_ord1, display = "bp")

# Agregamos las etiquetas y creamos una tabla
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Definimos la coordenadas para agregar las flechas
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

# Definimos las etiquetas para agregar las flechas
label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Creamos una nueva grafica
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




### Saving session info


data="a"
data2="B"

save.image("Toda_la_session.RData")
save(data, file="un_objeto.Rdata")


load("Toda_la_session.RData")
load("un_objeto.Rdata")


