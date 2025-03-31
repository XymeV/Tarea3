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
system.time(rarecurve(tabla, step = 100, sample = raremax, col = "red", label = FALSE ))

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
familia<- (tax_table(gb))
#menos de 5 lecturas en al menos 20% de muestras
gb
filter<- phyloseq::genefilter_sample(gb, filterfun_sample(function(x) x >= 5),
                                     A = 0.2*nsamples(gb))
gb_filter <- prune_taxa(filter, gb)
gb_filter

#Aglomerar a nivel de familia
aglomerado<-tax_glom(gb_filter, taxrank = "Family")
aglomerado

#Transformar a abundancias relativas
#Se calcuraron las abundancais relativas a nivel de familia, en un data frame que las muestra juntas
rabundancia<- tax_glom(GlobalPatterns, "Family") %>% transform_sample_counts(function(x) {
  x *100 / sum(x)
})
familia<-  data.frame(Family = tax_table(rabundancia)[,"Family"], Mean = rowMeans(otu_table(rabundancia)), row.names = NULL)
familia<- familia[order(-familia$Mean),]
familia

#Subset: solo Soil, Feces y Skin
sampledata<-sample_data(gb)
subset<- sampledata[sampledata$SampleType %in% c("Feces", "Soil", "Skin"), ]
gb_subset <- prune_samples(sample_names(gb) %in% rownames(subset), gb) 
View(subset)

#Diversidad alfa: Shannon, Simpson y Observed
gpalpha<- estimate_richness(gb_filter, measures = c("Shannon", "Simpson", "Observed"))
gpalpha
write.csv(gpalpha)
alfa<- plot_richness(gb, measures = c("Shannon", "Simpson", "Observed"))

#Boxplot de los tipos de muestra
tablagp<- as.data.frame(sample_data(gb)) #hacemos un data frame de estos datos, que son las muestras
tablas<- cbind(gpalpha, tablagp) #unimos los datos de los indices y las muestras

#Box plot  de la abundancia observada
ggplot(tablas, aes(x = SampleType, y =Observed, 
                   fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Diversidad (Observada)",
       x= "muestras",
       y = "abundancia")

#Box plot  de la abundancia Simpson
ggplot(tablas, aes(x = SampleType, y =Simpson, 
                   fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Diversidad  (Simpson)",
       x= "muestras",
       y = "abundancia")

#Box plot  de la abundancia Shannon
ggplot(tablas, aes(x = SampleType, y =Shannon, 
                   fill = SampleType)) +
  geom_boxplot() +
  labs(title = "Diversidad (Shannon)",
       x= "muestras",
       y = "abundancia")

#Prueba de Kruskal-Wallis
muestras <- sample_data(gb)$SampleType
tablas <- cbind(gpalpha, SampleType = muestras)

#Shannon
kruskals <- kruskal.test(Shannon ~ SampleType, data = tablas)
print(kruskals)

#Simpson
kruskalsimp <- kruskal.test(Simpson ~ SampleType, data = tablas)
print(kruskalsimp)

kruskalo <- kruskal.test(Observed ~ SampleType, data = tablas)
print(kruskalo)

#rank-Abundance
#volvemos a mencionar el subset anteriormente creado
sampledata <- sample_data(gb)
subset <- sampledata[sampledata$SampleType %in% c("Feces", "Soil", "Skin"), ]
gb_subset <- prune_samples(sample_names(gb) %in% rownames(subset), gb) 

#agrupar por familia
agb <- tax_glom(gb_subset, taxrank = "Family")

#transformar a abundancia relativa
rabundancia <- transform_sample_counts(agb, function(x) x * 100 / sum(x))


#otu table y taxa table
familytable <- otu_table(rabundancia)
familytabletaxa <- tax_table(rabundancia)

familia <- as.data.frame(as.table(familytable))
familia$SampleType <- rep(sample_data(rabundancia)$SampleType, each = nrow(tax_table(rabundancia)))
familia <- familia %>% filter(Freq > 0)
names(familia)[1] <- "Family"


#ordenar familias
familia <- familia %>%
  group_by(Family) %>%  
  mutate(TotalAbundance = sum(Freq)) %>%
  ungroup() %>%
  arrange(desc(TotalAbundance))

#familias por tipo de muestras
familia$Rank <- ave(familia$TotalAbundance, familia$SampleType, FUN = function(x) order(-x))

#grafica
ggplot(familia, aes(x = Rank, y = Freq)) +
  geom_line(aes(color = SampleType), size = 1) +  
  geom_point(aes(color = SampleType), size = 2) +  
  scale_x_continuous(trans = "log10") +  
  scale_y_continuous(trans = "log10") +  
  labs(title = "Rank-Abundance (Tipo de Muestra)",
       x = "Rango de Familia ",
       y = "Abundancia Relativa") +
  theme_minimal() +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    axis.text.y = element_text(size = 10)  
  ) +
  facet_wrap(~ SampleType, scales = "free_y") 


#Grafico apilado a nivel de filo
phylum<- tax_glom(gb, taxrank = "Phylum")
rabundancia<- transform_sample_counts(phylum, function(x) x * 100 / sum(x))

tablafilo<- otu_table(rabundancia)
taxafilo<- tax_table(rabundancia)

filo <- as.data.frame(as.table(tablafilo))

colnames(filo) <- c("SampleID", "Phylum", "Abundance")
filo<- filo %>% filter(Abundance > 0)

samplefilo <- as.data.frame(sample_data(rabundancia))
filo$SampleType <- samplefilo$SampleType[match(filo$SampleID, rownames(samplefilo))]


#grafico apilado
ggplot(filo, aes(x = SampleID, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(title = "Gráfico Apilado (Phylum)",
       x = "Muestras",
       y = "Abundancia Relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  facet_wrap(~ SampleType, scales = "free_x")

#con los 5 filos mas abundantes

todos <- filo %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(Abundance)) %>%
  arrange(desc(TotalAbundance))

cinco <- head(todos$Phylum, 5)
filtrar <- filo %>% filter(Phylum %in% cinco)

#Gráfico
ggplot(filtrar, aes(x = SampleID, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(title = "Gráfico Apilado de los 5 filos",
       x = "Muestras",
       y = "Abundancia Relativa") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotar las etiquetas del eje X
  facet_wrap(~ SampleType, scales = "free_x")

#Bray-Curtis distancia
otugb<- otu_table(gb)
otugb<- as(otugb, "matrix")
distanciabray<- vegdist(otugb, method = "bray")
distanciabray

#graficar
nmds <- metaMDS(distanciabray, k = 2, trymax = 100)
nmds_points <- as.data.frame(nmds$points)
nmds_points$SampleID <- rownames(nmds_points)

sample_data_df <- as.data.frame(sample_data(gb))
nmds_points$SampleType <- sample_data_df$SampleType[match(nmds_points$SampleID, rownames(sample_data_df))]

sum(is.na(nmds_points$SampleType))
head(nmds_points$SampleID) 
head(rownames(sample_data(gb)))

setdiff(nmds_points$SampleID, rownames(sample_data(gb)))
setdiff(rownames(sample_data(gb)), nmds_points$SampleID) 

nmds_points$SampleType <- sample_data(gb)$SampleType[match(nmds_points$SampleID, rownames(sample_data(gb)))]
sum(is.na(nmds_points$SampleType)) 

table(nmds_points$SampleID %in% rownames(sample_data(gb)))
issing_in_nmds <- setdiff(rownames(sample_data(gb)), nmds_points$SampleID)
missing_in_sample_data <- setdiff(nmds_points$SampleID, rownames(sample_data(gb)))

nmds_points$SampleID <- gsub("\\s+", "", nmds_points$SampleID)
nmds_points$SampleType <- sample_data(gb)$SampleType[match(nmds_points$SampleID, rownames(sample_data(gb)))]
sum(is.na(nmds_points$SampleType))

nmds_points$SampleType <- factor(nmds_points$SampleType)
ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = SampleType)) +
  geom_point(size = 3) +
  labs(title = "NMDS de Bray-Curtis para Diversidad Beta",
       x = "NMDS1", y = "NMDS2") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green", "purple"))

#Trate de corregir el codigo porque unicamente me aorecen NAs, revise y decia que podia ser porque no estaban correctamente asignados pero no pude cambiarlo.


#Permanova
sample_data_df <- as.data.frame(sample_data(gb))

rownames(sample_data_df) == rownames(otu_data) 
sample_data_df$SampleType <- factor(sample_data_df$SampleType)
permanova<- adonis(distanciabray ~ SampleType, sample_data_df, permutations = 999)
permanova


