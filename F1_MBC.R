#############################################################
### F1 Bioinformatics 
### Practical
### Part: Metabarcoding
### Lecturer: PD Dr. Alexander Keller
#############################################################
# This is material of a guided lecture. 
# Individual commands are explained while progressing through the script
# outputs are interpreted together on occurrence
# as well as intermediate files checked with "less" or "head"
#############################################################

#############################################################
# load packages
library(phyloseq)
library(ggplot2)
library(vegan)

# this package is not installed on the server globally, install locally then:
install.packages("bipartite", dependencies=T, repos="https://cloud.r-project.org")
library(bipartite)


#############################################################
# importing the data we just created
data.otu = otu_table(read.table("./combined.otu", sep="\t", header=T, row.names=1), taxa_are_rows=T)
data.tax = tax_table(as.matrix(read.table("./combined.tax", sep=",", header=T,row.names=1, fill=T)))
data.samp = import_qiime_sample_data("MapFile.txt")

# combining to a single object
dataset.comp = merge_phyloseq(data.otu, data.tax, data.samp)

#############################################################
# plotting diversity and species richness
pdf(file="plot_diversity.pdf")
	plot_richness(dataset.comp, x="BeeSpecies", measures=c("Shannon","Observed"))+geom_boxplot()
dev.off()

#############################################################
# transforming to relative data and removing species below 1% in samples
dataset.comp.rel.filter =transform_sample_counts(dataset.comp, function(x) x/sum(x))
otu_table(dataset.comp)[otu_table(dataset.comp.rel.filter)<0.01]<-0
otu_table(dataset.comp.rel.filter)[otu_table(dataset.comp.rel.filter)<0.01]<-0 
dataset.comp.rel.filter      = prune_taxa(taxa_sums(dataset.comp.rel.filter)>0, dataset.comp.rel.filter) # taxa die nach dem filtern gar nicht mehr vorkommen rausschmeißen

#############################################################
# ordination and plotting it
ordi = ordinate(dataset.comp.rel.filter, method="NMDS","bray", k=2)

pdf(file="plot_ordination.pdf")
	plot_ordination(dataset.comp.rel.filter, ordi,color="BeeSpecies")+geom_point(size=6)+theme_bw()
dev.off()


#############################################################
# plot bars and see species preferences
pdf(file="plot_barplot.pdf")
	plot_bar(dataset.comp.rel.filter,x="Family", fill="BeeSpecies")
dev.off()

#############################################################
# interaction network

dataset.comp.family<- tax_glom(dataset.comp.rel.filter,taxrank="Family")
taxa_names(dataset.comp.family) <- tax_table(dataset.comp.family)[,"Family"]

dataset.bees= merge_samples(dataset.comp.family, "BeeSpecies", fun=mean)

sample_names(dataset.bees)

pdf(file="plot_network.pdf")
	plotweb(round(t(data.frame(otu_table(dataset.bees)))*100))
dev.off()

