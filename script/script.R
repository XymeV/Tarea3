#script phyloseq

#PRIMERA PARTE

#Se cargan las librerias
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(dplyr)

#se carga un dataset de la libreria microbiome -> dietswap
data("dietswap", package = "microbiome")
ps<- dietswap
ps

#Conocer cualntas y cuales son las variables disponibles en los metadatos
head(sample_data(ps))


#calcular la curva de rarefaccion a partir del objeto phyloseq
#(Para esta parte me ayude de una publicación en Github)
data("dietswap")
taxa_are_rows(dietswap)
tabla<- t(otu_table(dietswap))
class(tabla)<- "matrix"
class(tabla)
tabla<- as(t(otu_table(dietswap)), "matrix")
class(tabla)
raremax<- min(rowSums(tabla))
system.time(rarecurve(tabla, step = 100, sample = raremax, col = "pink", label = FALSE ))

#Conocer la diversidad alfa del sitio mediante plot_richness.
#observed
plot_richness(ps, measures = "Observed")

#Shannon
plot_richness(ps, measures = "Shannon")

#Simpson
plot_richness(ps, measures = "Simpson")

  
#Filtrado y transformación
#Solo selecciionar los generos mas abundantes.
filtrado <- filter_taxa(ps, function(x) mean (x) >0.001, TRUE)
filtrado

#Diversidad Beta
ps1<-otu_table(ps)
ordencion<- ordinate(ps1, method="NMDS", distance="bray", formula=NULL)
ordencion 
#Se ocupó la tecnica de NMSDS tomada de: https://vaulot.github.io/tutorials/Phyloseq_tutorial.html#content
#la cual se usa cuando se tiene una matriz de disimilaridad o de distnacias.

#Gráfico PCoA
plot_ordination(ps, ordencion)

#Graficas apiladas por abundancia por taxon phyloseq
plot_bar(ps, fill = "Family") + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")
#Se agrupan por familia

#Curva Rank Abundance: para relaizar el cogido, me apoyé del puesto en https://www.researchgate.net/post/Plotting-rank-abundance-curve-in-phyloseq
#por el usuario Muhammad Arslan.

#Se transforma los taxones a porcentaje 
phyloTemp = transform_sample_counts(ps, function(x) 1e+02 * x/sum(x))

clusterData = psmelt(phyloTemp)
clusterData = dplyr::filter(clusterData, Abundance > 0)

#Se calcula la media y se escogen los taxones que se van a mostrar
clusterAgg = aggregate(Abundance ~ OTU + Family, data=clusterData, mean)
clusterAgg = aggregate(Abundance ~ OTU + Family,data=clusterData,mean)

#Filtrado y eleccion de los numeros de la tabla
clusterAgg = clusterAgg[order(-clusterAgg$Abundance),][1:100,]

ggplot(clusterAgg,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point(aes(color=Family), size = 3)+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  scale_y_log10()

#Tabla CSV
ps
csv<- as.data.frame(abundances(ps))
csv
write.csv(csv, "abundnacias.csv")


#SEGUNDA PARTE

#cargar los datos
data("GlobalPatterns")
gb<- GlobalPatterns
gb1<-otu_table(gb)

#menos de 5 lecturas en al menos 20% de muestras
gb
filter<- phyloseq::genefilter_sample(gb, filterfun_sample(function(x) x >= 5),
                                     A = 0.2*nsamples(gb))
gb_filter <- prune_taxa(filter, gb)
gb_filter

#Aglomerar a nivel de familia
gb = tax_glom(gb, taxrank = "Family")
gb

#Transformar a abundancias relativas


#Filtrar taxa con 

