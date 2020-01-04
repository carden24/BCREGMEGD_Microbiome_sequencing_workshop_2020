# Procesamiento de datos de microbioma usando *Mothur*
### Por Erick Cárdenas Poiré


---

# Tabla de contenidos
1. [Introducción](#p1)
2. [Preparativos](#p2)
    1. [Instalación de *Mothur*](#p2.1)
    2. [Archivos de referencia para *Mothur*](#p2.2)
    3. [Secuencias y archivos relacionados](#p2.3)
3. [Procesamiento inicial](#p3)
   1. [Concatenación de pares](#p3.1)
   2. [Remoción de secuencias anormales](#p3.2)
   3. [Dereplicación](#p3.3)
4. [Alineamiento](#p4)
5. [Eliminación de errores](#p5)
   1. [Reducción de ruido](#p5.1)
   2. [Remoción de quimeras](#p5.2)
   3. [Remoción de lineajes indeseados](#p5.3)
6. [Creación de tablas de OTUs](#p6)
7. [Creación de arboles filogenéticos](#p7)
8. [Exportación de datos](#p8)

## Introducción  <a name="p1"></a>
Este tutorial usa el programa *Mothur* para procesar archivos de secuenciación del gen de ARNr 16S. Este tutorial esta basado en el protocolo operativo estándar desarrollado por Patrick Schloss que se encuentra disponible [acá](http://www.mothur.org/wiki/MiSeq_SOP). 

**Requisitos:**
- Mothur instalado (ver abajo)
- Archivos de referencia de Mothur (en la carpeta MiSeq_SOP_files)
- Secuencias y archivos relacionados (en la carpeta MiSeq_SOP_files)

## Preparativos <a name="p2"></a>

Lo primero que hay que hacer es decidir donde poner todos los archivos para el tutorial. Lo ideal es crear un directorio de fácil acceso. Por ejemple podemos crear un directorio llamado *mothur* en el disco duro C.  

### Instalación de *Mothur* <a name="p2.1"></a>

La versión mas reciente de *Mothur* se encuentra [acá](https://github.com/mothur/mothur/releases/latest).  Para ejecutar este tutorial es necesario bajar el archivo y descomprimirlo dentro de la carpeta de trabajo. Esta versión de *Mothur* no tiene interfaz gráfica sino que se ejecuta tipeando comandos en el terminal de mothur. Luego de descomprimir el archivo haz click en el y se vera lo siguiente:

![Imagen de Mothur prompt](https://carden24.github.com/images/Mothur.jpg) 

Si es que esto no funciona se puede acceder a la terminal de windows directamente e iniciar *Mothur* desde ahí.

Para acceder al terminal desde windows apretar   **[Windows] + R**  
Luego tipear **cmd** y apretar **Enter**

Por predeterminado nos encontramos en el directorio del usuario actual ````c:\Users\erick````

Para movernos de directorio usamos el comando **cd**

````
cd ..             Cambiamos al directorio de arriba
cd c:\mothur      Cambiamos al directorio "mothur" en c:/

````

Cuando llegamos al directorio donde se encuentran los ejecutables de *Mothur* podemos iniciar el programa tipeando:
````
mothur.exe
````

En esta nueva terminal es que tenemos que tipear los comandos. *Mothur* mantiene una lista de los ultimos archivos creados asi que no es necesario muchas veces retipear el nombre del archivo sino solo reemplazar su nombre con "current". Para saber cuales son los archivos o parametros  que *Mothur* considera como los mas actuales se puede usar el comando *get.current()*. Para cambiar estos parametros se puede usar el comando *set.current*.

*Mothur* puede también recibir comandos directamente desde el terminal de Windows o ejecutar una serie de comandos escritos en un archivo de texto (un comando por linea). 


````
mothur "#fastq.info(fastq=test.fastq);get.current()"

mothur stability.batch
````

El archivo *stability.batch* contiene todos los comandos necesarios para el procesamiento estándar y se puede reutilizar para analizar otras muestras siempre y cuando se cambie el contenido del archivo *stability.files* que tiene la lista de archivos de secuencias.

### Archivos de referencia para *Mothur* <a name="p2.2"></a>

*Mothur* necesita archivos extras para poder funcionar. Los archivos se encuentran en la carpeta MiSeq_SOP_files. Es necesario bajarlos y descomprimirlos en la carpeta de trabajo. Estos son los archivos que vamos a usar para el tutorial:

- Alineamiento de referencia de *Silva* versión 102(Silva.bacteria.zip). Útil para tutorial, versiones mas recientes disponibles [acá](https://www.mothur.org/wiki/Silva_reference_files).
- Taxonomía del RDP formateado para *Mothur* (Trainset14_032015.pds.zip). Mothur mantiene la lista de taxonomías del RDP [acá](https://www.mothur.org/wiki/RDP_reference_files). *Mothur recomienda* usar la taxonomia de Greengenes. *Mothur* mantiene una lista de las  versiones mas recientes [acá](https://www.mothur.org/wiki/Greengenes-formatted_databases).


### Secuencias y archivos relacionados <a name="p2.3"></a>

Finalmente necesitamos bajar las secuencias de Miseq y archivos relacionados al experimento. El tutorial original usa 21 pares, esta usa 9 pares para hacer el tutorial mas corto.

- Secuencias en formato fastq (Miseq SOP.zip)
- Diseño experimental (mouse.time.design)
- Metadatos (mouse.dpw.metadata)
- stability.files (lista de archivos de secuencias)

## Procesamiento inicial <a name="p3"></a>

### Concatenación de pares <a name="p3.1"></a>

El comando *make.contigs* une el par de secuencias que provienen de la misma molécula. Los primers que se usan para amplificar la region V4 son los 515F y 806R por lo que el producto deber ser de unos 290 bases y al secuenciar este producto por ambos lados con secuencias de 250 bases crea un gran solapamiento. *Mothur* une las secuencias y usa los valores de calidad para asignar bases en la región que se sobrelapa.  


```
make.contigs(file=stability.files, processors=2)
```

Este comando también genera archivos que se necesitan después. 
:
  
- *stability.trim.contigs.fasta* : Nuevas secuencias unidas
- *stability.contigs.groups* : El grupo al cual pertenece cada secuencias
- *stability.contigs.report* : Reporte del proceso de concatenamiento
  
Para crear estadisdicas de las secuencias usamos el comando *summary.seqs*


```
summary.seqs()
 
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       248     248     0       3       1
2.5%-tile:      1       252     252     0       3       3810
25%-tile:       1       252     252     0       4       38091
Median:         1       252     252     0       4       76181
75%-tile:       1       253     253     0       5       114271
97.5%-tile:     1       253     253     6       6       148552
Maximum:        1       502     502     249     243     152360
Mean:           1       252.811 252.811 0.70063         4.44854
# of Seqs:      152360
 ```
El resultado muestra que la mayoría de secuencias tienen entre 248 y 253 bases, y hay algunas secuencias que tienen hasta 502 bases lo cual sugiere que probablemente no se unieron bien ya que esperamos un producto de 290 bases como maximo.  

 
### Remoción de secuencias anormales <a name="p3.2"></a>

Usamos el comando *screen.seqs*  para remover secuencias de acuerdo a su tamaño y el numero de bases ambiguas en ellas (Ns). Este paso remueve secuencias erróneas y artefactos. En este caso removemos secuencias con bases ambiguas y secuencias mas largas que 275 bp. Este numero depende del tamaño del la región definida por los primers.


 
```
screen.seqs(fasta=current, group=current, summary=current, maxambig=0, maxlength=275)


>summary.seqs()
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       250     250     0       3       1
2.5%-tile:      1       252     252     0       3       3222
25%-tile:       1       252     252     0       4       32219
Median:         1       252     252     0       4       64437
75%-tile:       1       253     253     0       5       96655
97.5%-tile:     1       253     253     0       6       125651
Maximum:        1       270     270     0       12      128872
Mean:   1       252.462 252.462 0       4.36693
# of Seqs:      128872
```
 

### Dereplicación <a name="p3.3"></a>

En este paso removemos secuencias idénticas para reducir la carga en la computadora. *Mothur* se encarga de acordarse de que hizo esto. 
Este es un paso no cambia la calidad de las secuencias pero reduce la carga computacional (tiempo de procesamiento y memoria requerida).

```
unique.seqs(fasta=current)
``` 

El protocolo requiere que usemos *count.seqs* para crear una tabla que registra las secuencias repetidas.

```
count.seqs(name=current, group=current)

```

## Alineamiento <a name="p4"></a> 

En este paso alineamos nuestras secuencias con el alineamiento de referencia de SILVA, un alineamiento curado manualmente de gran calidad. Este es un archivo grande con casi 15000 secuencias y mas de 50000 posiciones. Para hacer nuestra trabajo mas fácil vamos a editar el alineamiento de Silva a las región que nos interesa (V4) y luego alinearemos nuestras secuencias con esta selección. Finalmente editaremos el alineamiento que incluye nuestras secuencias.  


El primer comando, *pcr.seqs*, edita el alineamiento de Silva a nuestra región de interés. Si es que no sabes las coordinadas de la región de interés, saltea el tutorial hasta el comando *align.seqs* y ahí usa la opción *reference=silva.bacteria.fasta*.

```
pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F)
```

Ahora podemos renombrar el nuevo archivo creado con algo mas fácil de entender usando el comando *rename.file*. 

```
rename.file(input=silva.bacteria.pcr.fasta, new=silva.v4.fasta)
```

Ahora podemos alinear las secuencias al alineamiento maestro de Silva para la region V4. 


```
align.seqs(fasta=current, reference=silva.v4.fasta)
```

Luego de alinear las secuencias podemos ejecutar *summary.seqs()* de nuevo  para ver estadísticas sobre distribución de las secuencias en el alineamiento. 

````
summary.seqs(count=current)

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1250    10693   250     0       3       1
2.5%-tile:      1968    11550   252     0       3       3222
25%-tile:       1968    11550   252     0       4       32219
Median:         1968    11550   252     0       4       64437
75%-tile:       1968    11550   253     0       5       96655
97.5%-tile:     1968    11550   253     0       6       125651
Maximum:        1982    13400   270     0       12      128872
Mean:           1967.99 11550   252.462 0       4.36693
# of unique seqs:       16426
total # of seqs:        128872
````

Al evaluar el alineamiento podemos ver que la mayoría de las secuencias empieza en la posición 1968 y termina en la posición 11550. Esta evaluación es útil detectar secuencias con muchas inserciones o que empiezan o terminan muy tarde. En este paso también podemos remover secuencias con muchos homopolímeros ya que estas tienden a ser erróneas.  

Utilizamos entonces el comando *screen.seqs* para eliminar secuencias que empiezan mucho antes o terminan mucho despues ademas de secuencias con muchos homopolímeros.
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

Ahora podemos editar el alineamiento y  remover posiciones que solo tienen gaps ya que no contribuyen con datos. El comando *screen.seqs* también necesita saber cual es el carácter que el alineamiento usa para indicar que no hay datos (trump character). En nuestro caso el alineamiento de Silva usa ".".

```
filter.seqs(fasta=current, vertical=T, trump=.)
```

Ejecutamos *unique.seqs* una vez mas ya que algunas secuencias pueden ser idénticas luego del alineamiento. 

```
unique.seqs(fasta=current, count=current)
```

Finalmente ejecutamos summary.seqs() y vemos que el alineamiento tiene muchas menos posiciones lo cual hace el calculo posterior mas fácil. 

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

## Eliminación de errores <a name="p5"></a>

Estos pasos reducen mejoran la calidad de los datos al eliminar errores de secuenciación, secuencias provenientes del anfitrión, y algunos artefactos generados por la PCR.  

### Reducción de ruido <a name="p5.1"></a>

El primer paso es agrupar secuencias altamente similares con *pre.cluster*. Este comando primero agrupa secuencias que difieran en unas cuantas posiciones, este caso no mas de dos diferencias (1 por cada 100 bases). Si hay varias secuencias similar el algoritmo escoge a la secuencia mas abundante como la secuencias correcta y asume que las demás son diferentes debido a errores de secuenciación. Este paso no remueve secuencias.

```
pre.cluster(fasta=current, count=current, diffs=2)
```

### Remoción de quimeras <a name="p5.2"></a>

Las quimeras son artefactos de la PCR creados cuando la secuencias de una especie se alinean con las de otra especies. El resultado de este encuentro casi amoroso es un artefacto que no representa la diversidad real de la comunidad.

El primer paso detecta las quimeras. El segundo las remueve. Deben ejecutarse uno después del otro   
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
vsearch removió el 8.2% de las secuencias (n=10511), y se redujo además el numero de secuencias únicas lo cual hace mas fácil el proceso.  


### Remoción de lineajes indeseados <a name="p5.3"></a>

En algunos casos nuestra PCR puede amplificar secuencias de organelos del huésped (mitocondrias y cloroplastos), secuencias de arqueas o eucariotas (con menos especificidad y poca sensibilidad) y otros productos no especificos. Estas secuencias deben ser detectadas y removidas pues no representan a la comunidad microbiana.  
El primer paso es clasificar todas las secuencias de la comunidad con *classify.sequences* y luego remover los lineajes no deseados con *remove.lineage*. El primer paso usa el método Bayesiano de clasificación, y la taxonomica de referencia del RDP. Esta taxonomía es útil para el tutorial pero se recomienda que para trabajos reales se use la taxonomía de *Greengenes* que clasifica hasta nivel de especies disponible [acá](https://www.mothur.org/wiki/Greengenes-formatted_databases). 


```
classify.seqs(fasta=current, count=current, reference=trainset14_032015.pds.fasta, taxonomy=trainset14_032015.pds.tax, cutoff=80)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-mitochondria-unknown-Archaea-Eukaryota)


Si es que trabajamos con Greengenes
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

El tutorial de Miseq original explica como utilizar comunidades modelos (con abundancias definidas) para calcular tasas de error. Si es que planeas usar este tipo de control positivo puedes ver con mas detalles el análisis en el [tutorial original de Mothur](https://www.mothur.org/wiki/MiSeq_SOP). En este caso simplemente vamos a remover las  muestras que de comunidades modelos utilizando el comando *remove.seqs*. 


```
remove.groups(count=current, fasta=current, taxonomy=current, groups=Mock)
```

## Creación de tablas de OTUS <a name="p6"></a>

El siguiente objetivo es agrupar las secuencias en unidades taxonómicas operacionales (OTUs), grupos de secuencias definidos por la similitud entre ellas. La forma tradicional es primero crear una matriz de distancias con *dist.seqs* y luego crear grupos de secuencias con esta matriz con *cluster*. Este método es muy lento cuando se trabajan con miles de secuencias ya que la matriz ocupa mucho espacio en memoria. 

````
No ejecutar estos comandos !
dist.seqs(fasta=current, cutoff=0.03)
cluster(column=current, count=current)

````

La alternativa practica es usar *cluster.split*. Este método primero crea grupos basándose en la matriz de distancia o en la clasificación taxonómica de las secuencias y posteriormente aplica los métodos de agrupamientos sobre las selecciones individuales. Es este caso usaremos la clasificación taxonómica para separar los grupos a nivel de orden (nivel taxonómico), y luego aplicaremos el método algoritmo *opticlust* para crear OTUs. El comando *cluster.split* necesita una valor para saber hasta que porcentaje de similitud hay que crear las matrices de distancia. Para este algoritmo, cutoff=0.03 es suficiente, para el método *average neighbour* se recomiendan valores mas altos que el parametro final (usar 0.15 si es que se quiere trabajar a nivel de especies ).


```
cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.03)
```

Ahora podemos usar las matrices de distancia para crear las tablas de OTUs. En este caso nos interesa trabajar a nivel de especie que es aproximadamente 3%.


```
make.shared(list=current, count=current, label=0.03)
```
La tabla que se genero es lo necesario para hacer los análisis de diversidad. 

El ultimo paso del protocolo es asignar una clasificación taxonómica a cada OTU. Para esto usaremos el comando *classify.otu* que genera para la clasificación de consenso para el OTU. 

```
classify.otu(list=current, count=current, taxonomy=current, label=0.03)
```

## Creación de arboles filogenético <a name="p7"></a>

Algunos métodos como Unifrac requieren saber la localización de cada OTU en un árbol filogenético de este proyecto. Este proceso se basa en el alineamiento y crea primero una matriz de distancias y luego usa esa matriz para crear un árbol basado en distancias (método Neighbour-joining).

```
dist.seqs(fasta=current, output=lt)
clearcut(phylip=current)
```

## Exportación de datos <a name="p8"></a>

Mothur genera muchos archivos intermedios con nombres complicados. La forma mas fácil de transferir estos archivos a otros programas es convirtiéndolo a un objeto en formato biom. Este formato es una forma eficiente y facil de procesar (para las computadoras). Para mas información de este formato ver [acá] http://biom-format.org/documentation/biom_format.html. 


````
make.biom(shared=current, constaxonomy=current, metadata=mouse.dpw.metadata)

````  


## Análisis de filotipos   [Opcional]

Si es que no nos interesa crear tablas de OTUs, especies definidas por similitud, sino trabajar con tablas de especies definidas por su clasificación taxonómica (filotipos) podemos usar los siguientes tres comandos con los datos ya procesados (sin artefactos y ruido).

````
phylotype(taxonomy=current, label=1);
make.shared(list=current, count=current);
classify.otu(list=current, taxonomy=current);
````

El primer comando agrupa las secuencias según su taxonomía, la opción "label=1" especifica que vamos a trabajar con la taxonomía mas especifica. Si es que no especificamos este parámetro el resultado va crear listas de secuencias a diferente niveles de taxonomía. 
El segundo comando crea la tabla de filotipos x muestras. El último comando crea otra tabla con la taxonomía de consenso del filotipo.
Hay que ser cuidadoso en escoger la tabla con la cual trabajar ya que tienen nombres similares.

```
Archivos "Shared":
stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.shared  <-Basado en OTUs
stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.shared           <-Basado en filotipos

Archivos de taxonomia de consenso
stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.1.cons.tax.summary              <-Basado en OTUs
stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.tax.summary  <-Basado en filotipos

```

**Recomendación final**
Una vez terminado el procesamiento hay que grabar los fastq originales, la lista de comando que usamos, los logs del programa (que tienen los comandos y los resultados), las tablas de OTUs (archivos SHARED), la clasificación de los OTUs (.cons.tax.summary), los metadatos, y el árbol filogenetico (archivo .phylip.tre) si es que nos interesa usar indices como UNIFRAC. Todos los demás resultados se puede borrar o comprimir.
 










