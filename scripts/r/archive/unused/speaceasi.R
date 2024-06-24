library(remotes)
#devtools::install_github("zdk123/SpiecEasi") #had to uninstall rcpp and install it again (update not)
library(SpiecEasi)
library(phyloseq)
library(metagMisc)
library(igraph)

getwd()
library(SpiecEasi)

##################general##############

## Imported otu table 
Seed_Bac_otu <- read.table("tableBAC_130-150_original.txt", header=TRUE, row.names="OTUID"); Seed_Bac_otu
#Seed_Bac_otu <- t(Seed_Bac_otu)
Seed_Bacteria_otu <- otu_table(Seed_Bac_otu, taxa_are_rows = TRUE); Seed_Bacteria_otu
## Imported taxonomy
Seed_bacteria_tax <- read.table("taxonomy.txt", header=TRUE, row.names="OTUID")
Seed_bacteria_tax <- as.matrix(Seed_bacteria_tax)
Seed_bacteria_tax <- tax_table(Seed_bacteria_tax)
## Imported metadata file
seed_metadata <- read.csv("Seed_metadata.csv")
rownames(seed_metadata) <- seed_metadata[,1]
seed_metadata <- sample_data(seed_metadata); seed_metadata

#Convert to phyloseq object
Bacteria_seed <- merge_phyloseq(Seed_Bacteria_otu,Seed_bacteria_tax, seed_metadata);Bacteria_seed

#Subset taxa if needed 
Bacteria_seed = subset_taxa(Bacteria_seed, Kingdom=="Bacteria"); Bacteria_seed
#Bacteria_seed <- prune_taxa(taxa_sums(Bacteria_seed)>0, Bacteria_seed); Bacteria_seed

Bacteria_seed = subset_samples(Bacteria_seed, country != "Austria"); Bacteria_seed


Bacteria_seed_otu = as(otu_table(Bacteria_seed), "matrix")
# Coerce to data.frame
Bacteria_seed_otu = as.data.frame(Bacteria_seed_otu)
write.csv(Bacteria_seed_otu ,file = "Bacteria_seed_otu.csv")




sums_all<- sample_sums(Bacteria_seed); sums_all
sums_all <- data.frame(sums_all)
#sums_matrix_BACTERIA <- as.matrix(sums_BACTERIA)
colSums(sums_all) #25271878     
#6127107/27808166=0.22
max(sums_all) #565382
min(sums_all) #25
sums_all

#write.csv(sums_all_df ,file = "sums_all_df.csv")

sums_all%>%
  arrange(sums_all)


Bacteria_seed = subset_taxa(Bacteria_seed, Order!="Chloroplast"); Bacteria_seed 
Bacteria_seed = subset_taxa(Bacteria_seed, Family!="Mitochondria"); Bacteria_seed



sums_BACTERIA <- sample_sums(Bacteria_seed); sums_BACTERIA
sums_BACTERIA <- data.frame(sums_BACTERIA)
#sums_matrix_BACTERIA <- as.matrix(sums_BACTERIA)
colSums(sums_BACTERIA) #5253461   
#6127107/27808166=0.22
max(sums_BACTERIA) #435504
min(sums_BACTERIA) #22
sums_BACTERIA

sums_BACTERIA%>%
  arrange(sums_BACTERIA)

################################### Decontam #############

###run1
Seed1 = subset_samples(Bacteria_seed, Sequencing_Run=="Run1"); Seed1

sample_data(Seed1)$is.neg <- sample_data(Seed1)$Sample_or_Control == "Control"
contamdf.prev_Seed1 <- isContaminant(Seed1, method="prevalence", neg="is.neg")
table(contamdf.prev_Seed1$contaminant)

head(which(contamdf.prev_Seed1$contaminant))

contamdf.prev05_seed1 <- isContaminant(Seed1, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05_seed1$contaminant)


#ps.pa_seed1 <- transform_sample_counts(Seed1, function(abund) 1*(abund>0))
ps.pa.neg_seed1 <- prune_samples(sample_data(Seed1)$Sample_or_Control == "Control", Seed1)
ps.pa.pos_seed1 <- prune_samples(sample_data(Seed1)$Sample_or_Control == "Sample", Seed1)
df.pa_seed1 <- data.frame(pa.pos_seed1=taxa_sums(ps.pa.pos_seed1), pa.neg_seed1=taxa_sums(ps.pa.neg_seed1),
                          contaminant=contamdf.prev05_seed1$contaminant)
ggplot(data=df.pa_seed1, aes(x=pa.neg_seed1, y=pa.pos_seed1, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

ps_seed1_decontam_05 <- prune_taxa(!contamdf.prev05_seed1$contaminant, Seed1)
ps_seed1_decontam_05

ps_seed1_decontam_01 <- prune_taxa(!contamdf.prev_Seed1$contaminant, Seed1)
ps_seed1_decontam_01

###run2


Seed2 = subset_samples(Bacteria_seed, Sequencing_Run=="Run_H3"); Seed2
Seed2<- subset_samples(Seed2, Sample_Index != "ctrl_118")
Seed2<- subset_samples(Seed2, Sample_Index != "ctrl_120")
Seed2<- subset_samples(Seed2, Sample_Index != "ctrl_121") #ctrls associated with rye samples


sample_data(Seed2)$is.neg <- sample_data(Seed2)$Sample_or_Control == "Control"

contamdf.prev_Seed2 <- isContaminant(Seed2, method="prevalence", neg="is.neg")
table(contamdf.prev_Seed2$contaminant)
head(which(contamdf.prev_Seed2$contaminant))

contamdf.prev05_seed2 <- isContaminant(Seed2, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05_seed2$contaminant)
head(which(contamdf.prev05_seed2$contaminant))


ps.pa_seed2 <- transform_sample_counts(Seed2, function(abund) 1*(abund>0))
ps.pa.neg_seed2 <- prune_samples(sample_data(ps.pa_seed2)$Sample_or_Control == "Control", ps.pa_seed2)
ps.pa.pos_seed2 <- prune_samples(sample_data(ps.pa_seed2)$Sample_or_Control == "Sample", ps.pa_seed2)
df.pa_seed2 <- data.frame(pa.pos_seed2=taxa_sums(ps.pa.pos_seed2), pa.neg_seed2=taxa_sums(ps.pa.neg_seed2),
                          contaminant=contamdf.prev_Seed2$contaminant)
ggplot(data=df.pa_seed2, aes(x=pa.neg_seed2, y=pa.pos_seed2, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

ps_seed2_decontam_05 <- prune_taxa(!contamdf.prev05_seed2$contaminant, Seed2)
ps_seed2_decontam_05

ps_seed2_decontam_01 <- prune_taxa(!contamdf.prev_Seed2$contaminant, Seed2)
ps_seed2_decontam_01

#merge phyloseq objects
seed_decontaminated_05 <- merge_phyloseq(ps_seed1_decontam_05, ps_seed2_decontam_05);seed_decontaminated_05

#seed_decontaminated_01 <- merge_phyloseq(ps_seed1_decontam_01, ps_seed2_decontam_01)


seed_decontaminated_05<- subset_samples(seed_decontaminated_05, country != "ctrl")
seed_decontaminated_05<- subset_samples(seed_decontaminated_05, country != "Extrctrl")

seed_decontaminated_05<- subset_samples(seed_decontaminated_05, identity != "DFR21CE")

#seed_decontaminated_01<- subset_samples(seed_decontaminated_01, country != "ctrl")
#seed_decontaminated_01<- subset_samples(seed_decontaminated_01, country != "Extrctrl")


#seed_decontaminated_genus = tax_glom(seed_decontaminated_05, taxrank="Genus")

#seed_decontaminated_05_country_filtered <- filter_taxa(seed_decontaminated_05_country, function (x) {sum(x > 0) > 1}, prune=TRUE) #filter singeltons


sums_BACTERIA_decontam <- sample_sums(seed_decontaminated_05); sums_BACTERIA_decontam
sums_BACTERIA_decontam <- data.frame(sums_BACTERIA_decontam)
#sums_matrix_BACTERIA <- as.matrix(sums_BACTERIA)
colSums(sums_BACTERIA_decontam) #4522046   
#6127107/27808166=0.22
max(sums_BACTERIA_decontam) #434850
min(sums_BACTERIA_decontam) #22
sums_BACTERIA_decontam

write.csv(sums_BACTERIA_decontam ,file = "sums_BACTERIA_decontam.csv")

sums_BACTERIA_decontam%>%
  arrange(sums_BACTERIA_decontam)

############ Subset ####

seed_decontaminated_05_country = subset_samples(seed_decontaminated_05, funfact =="country"); seed_decontaminated_05_country 
seed_decontaminated_05_breeding = subset_samples(seed_decontaminated_05, funfact =="breeding"); seed_decontaminated_05_breeding 

#seed_decontaminated_05_country = subset_samples(seed_decontaminated_05_country, identity !="DFR21CE"); seed_decontaminated_05_country 

seed_decontaminated_05_ctrl = subset_samples(seed_decontaminated_05, Sample_or_Control =="Control"); seed_decontaminated_05_ctrl 


###################### all samples as network ################


Bacteria_seed_breeding_filtered <- phyloseq_filter_prevalence(seed_decontaminated_05_breeding, prev.trh = 0.01);Bacteria_seed_breeding_filtered  # 724
se.mb.breeding_all_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered, method='mb', lambda.min.ratio=1e-3, verbose=TRUE, nlambda=30, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_all_0.1_30) #0.04986547
sum(getRefit(se.mb.breeding_all_0.1_30))/2 # 179

betaMat_all=as.matrix(symBeta(getOptBeta(se.mb.breeding_all_0.1_100)))

positive_all=length(betaMat_all[betaMat_all>0.3])/2; positive_all #41
negative_all=length(betaMat_all[betaMat_all<-0.3])/2; negative_all #0

ig2.mb_breeding_all <- adj2igraph(getRefit(se.mb.breeding_all_0.1_100),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered)))
plot_network(ig2.mb_breeding_all, Bacteria_seed_breeding_filtered, type='taxa', color = "Phylum")

##### per Breeding cycle merged

breeding_C3 = subset_samples(seed_decontaminated_05_breeding, Breeding=="C3"); breeding_C3
breeding_filtered_C3 <- phyloseq_filter_prevalence(breeding_C3, prev.trh = 0.1, abund.trh = NULL);breeding_filtered_C3  # 64 taxa

se.mb.breeding_filtered_C3_default <- spiec.easi(breeding_filtered_C3, method='mb', verbose=TRUE, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #  Optimal lambda may be larger than the supplied values
se.mb.breeding_filtered_C3_0.01_20 <- spiec.easi(breeding_filtered_C3, method='mb',  lambda.min.ratio=0.01, nlambda= 20, verbose=TRUE, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_filtered_C3_0.01_20) #0.04986547
sum(getRefit(se.mb.breeding_filtered_C3_0.01_20))/2 # 38
se.mb.breeding_filtered_C3_0.01_20$select$stars$summary

betaMat_all=as.matrix(symBeta(getOptBeta(se.mb.breeding_filtered_C3_0.01_20)))

positive_all=length(betaMat_all[betaMat_all>0])/2; positive_all #41
negative_all=length(betaMat_all[betaMat_all<0])/2; negative_all #0

ig2.mb_breeding_C3 <- adj2igraph(getRefit(se.mb.breeding_filtered_C3_0.01_20),  vertex.attr=list(name=taxa_names(breeding_filtered_C3)))
plot_network(ig2.mb_breeding_C3, breeding_filtered_C3, type='taxa', color = "Phylum")



breeding_C5 = subset_samples(seed_decontaminated_05_breeding, Breeding=="C5"); breeding_C5
breeding_filtered_C5 <- phyloseq_filter_prevalence(breeding_C5, prev.trh = 0.1, abund.trh = NULL);breeding_filtered_C5  # 70 taxa

se.mb.breeding_filtered_C5_default <- spiec.easi(breeding_filtered_C5, method='mb', verbose=TRUE, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #  Optimal lambda may be larger than the supplied values
se.mb.breeding_filtered_C5_0.01_20 <- spiec.easi(breeding_filtered_C5, method='mb',  lambda.min.ratio=0.01, nlambda= 20, verbose=TRUE, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_filtered_C5_0.01_20) #0.04150261
sum(getRefit(se.mb.breeding_filtered_C5_0.01_20))/2 # 38
se.mb.breeding_filtered_C5_0.01_20$select$stars$summary

betaMat_all=as.matrix(symBeta(getOptBeta(se.mb.breeding_filtered_C3_0.01_20)))

positive_all=length(betaMat_all[betaMat_all>0])/2; positive_all #41
negative_all=length(betaMat_all[betaMat_all<0])/2; negative_all #0

ig2.mb_breeding_C3 <- adj2igraph(getRefit(se.mb.breeding_filtered_C3_0.01_20),  vertex.attr=list(name=taxa_names(breeding_filtered_C3)))
plot_network(ig2.mb_breeding_C3, breeding_filtered_C3, type='taxa', color = "Phylum")

breeding_C7 = subset_samples(seed_decontaminated_05_breeding, Breeding=="C7"); breeding_C7
breeding_filtered_C7 <- phyloseq_filter_prevalence(breeding_C7, prev.trh = 0.1, abund.trh = NULL);breeding_filtered_C7  # 70 taxa

se.mb.breeding_filtered_C5_default <- spiec.easi(breeding_filtered_C5, method='mb', verbose=TRUE, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #  Optimal lambda may be larger than the supplied values
se.mb.breeding_filtered_C7_0.01_20 <- spiec.easi(breeding_filtered_C7, method='mb',  lambda.min.ratio=0.01, nlambda= 20, verbose=TRUE, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_filtered_C7_0.01_20) #0.04150261
sum(getRefit(se.mb.breeding_filtered_C7_0.01_20))/2 # 38
se.mb.breeding_filtered_C7_0.01_20$select$stars$summary


breeding_C8 = subset_samples(seed_decontaminated_05_breeding, Breeding=="C8"); breeding_C8
breeding_filtered_C8 <- phyloseq_filter_prevalence(breeding_C8, prev.trh = 0.1, abund.trh = NULL);breeding_filtered_C8 # 70 taxa

se.mb.breeding_filtered_C5_default <- spiec.easi(breeding_filtered_C5, method='mb', verbose=TRUE, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #  Optimal lambda may be larger than the supplied values
se.mb.breeding_filtered_C8_0.01_20 <- spiec.easi(breeding_filtered_C8, method='mb',  lambda.min.ratio=0.01, nlambda= 20, verbose=TRUE, sel.criterion='stars', pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_filtered_C8_0.01_20) #0.04150261
sum(getRefit(se.mb.breeding_filtered_C8_0.01_20))/2 # 38
se.mb.breeding_filtered_C8_0.01_20$select$stars$summary

ig2.mb_breeding_C8 <- adj2igraph(getRefit(se.mb.breeding_filtered_C8_0.01_20),  vertex.attr=list(name=taxa_names(breeding_filtered_C8)))
plot_network(ig2.mb_breeding_C8, breeding_filtered_C8, type='taxa', color = "Phylum")

#########################

Bacteria_seed_breeding_C3 = subset_samples(seed_decontaminated_05_breeding, identity=="AC3"); Bacteria_seed_breeding_C3
Bacteria_seed_breeding_C5 = subset_samples(seed_decontaminated_05_breeding, identity=="AC5"); Bacteria_seed_breeding_C5
Bacteria_seed_breeding_C7 = subset_samples(seed_decontaminated_05_breeding, identity=="AC7"); Bacteria_seed_breeding_C7
Bacteria_seed_breeding_C8 = subset_samples(seed_decontaminated_05_breeding, identity=="AC8"); Bacteria_seed_breeding_C8

Bacteria_seed_breeding_C3_BE = subset_samples(seed_decontaminated_05_breeding, identity=="C3_BE"); Bacteria_seed_breeding_C3_BE
Bacteria_seed_breeding_C5_BE = subset_samples(seed_decontaminated_05_breeding, identity=="C5_BE"); Bacteria_seed_breeding_C5_BE
Bacteria_seed_breeding_C7_BE = subset_samples(seed_decontaminated_05_breeding, identity=="C7_BE"); Bacteria_seed_breeding_C7_BE
Bacteria_seed_breeding_C8_BE = subset_samples(seed_decontaminated_05_breeding, identity=="C8_BE"); Bacteria_seed_breeding_C8_BE

Bacteria_seed_breeding_filtered_C3 <- phyloseq_filter_prevalence(Bacteria_seed_breeding_C3, prev.trh = 0.1, abund.trh = NULL);Bacteria_seed_breeding_filtered_C3  # 185 taxa
Bacteria_seed_breeding_filtered_C5 <- phyloseq_filter_prevalence(Bacteria_seed_breeding_C5, prev.trh = 0.2, abund.trh = NULL);Bacteria_seed_breeding_filtered_C5  # 191 taxa
Bacteria_seed_breeding_filtered_C7 <- phyloseq_filter_prevalence(Bacteria_seed_breeding_C7, prev.trh = 0.2, abund.trh = NULL);Bacteria_seed_breeding_filtered_C7  # 94 taxa
Bacteria_seed_breeding_filtered_C8 <- phyloseq_filter_prevalence(Bacteria_seed_breeding_C8, prev.trh = 0.2, abund.trh = NULL);Bacteria_seed_breeding_filtered_C8  # 115 taxa

Bacteria_seed_breeding_filtered_C3_BE <- phyloseq_filter_prevalence(Bacteria_seed_breeding_C3_BE, prev.trh = 0.2, abund.trh = 0.001);Bacteria_seed_breeding_filtered_C3_BE  # 69 taxa
Bacteria_seed_breeding_filtered_C5_BE <- phyloseq_filter_prevalence(Bacteria_seed_breeding_C5_BE, prev.trh = 0.2, abund.trh = 0.01);Bacteria_seed_breeding_filtered_C5_BE  # 110 taxa
Bacteria_seed_breeding_filtered_C7_BE <- phyloseq_filter_prevalence(Bacteria_seed_breeding_C7_BE, prev.trh = 0.2, abund.trh = 0.01);Bacteria_seed_breeding_filtered_C7_BE  # 156 taxa
Bacteria_seed_breeding_filtered_C8_BE <- phyloseq_filter_prevalence(Bacteria_seed_breeding_C8_BE, prev.trh = 0.2, abund.trh = 0.01);Bacteria_seed_breeding_filtered_C8_BE  # 166 taxa

Bacteria_seed_breeding_all <- phyloseq_filter_prevalence(seed_decontaminated_05_breeding, prev.trh = 0.000001, abund.trh = NULL);Bacteria_seed_breeding_unfiltered_C3  # 185 taxa
se.mb.breeding_unf_C3_0.01_100 <- spiec.easi(Bacteria_seed_breeding_C3, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=100, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_unf_C3_0.1_100) #0.04900493
sum(getRefit(se.mb.breeding_unf_C3_0.1_100))/2 # 14
se.mb.breeding_unf_C3_0.1_100$select$stars$summary #0.01615685    


#Bacteria_seed_breeding_filtered_C3_BE_otu = as(otu_table(Bacteria_seed_breeding_filtered_C3_BE), "matrix")

bargs <- list(rep.num=50, seed=10010, conffile="parallel")

library(batchtools)

#se.mb.breeding_C3_0.01_100 <- spiec.easi(Bacteria_seed_breeding_filtered_C3, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=100, sel.criterion='stars',pulsar.select='batch', pulsar.params=bargs) #abort
#getStability(se.mb.breeding_C3_0.01_100) #0.04986547 #45
#sum(getRefit(se.mb.breeding_C3_0.01_100))/2 # 278
#se.mb.breeding_C3_0.02_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C3, method='mb', lambda.min.ratio=0.02, verbose=TRUE, nlambda=50, sel.criterion='stars',pulsar.select='batch', pulsar.params=bargs) #abort
#getStability(se.mb.breeding_C3_0.02_50) #0.04986547 #46 #smallbin double linked list corrupted
#sum(getRefit(se.mb.breeding_C3_0.02_50))/2 # 278 #46
#se.mb.breeding_C3_0.1_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C3, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=10, sel.criterion='stars',pulsar.select='batch', pulsar.params=bargs) #abort
#getStability(se.mb.breeding_C3_0.1_10) #0.04106905 
#sum(getRefit(se.mb.breeding_C3_0.1_10))/2 # 255 #44
#se.mb.breeding_C3_0.1_10$select$stars$summary
#se.mb.breeding_C3_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C3, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, sel.criterion='stars',pulsar.select='batch', pulsar.params=bargs) #abort
#getStability(se.mb.breeding_C3_0.1_30) #0.04399089 
#sum(getRefit(se.mb.breeding_C3_0.1_30))/2 # 230 #34
se.mb.breeding_C3_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C3, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, sel.criterion='stars',pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C3_0.1_50) #0.04794267 
sum(getRefit(se.mb.breeding_C3_0.1_50))/2 # 233 #34
se.mb.breeding_C3_default<- spiec.easi(Bacteria_seed_breeding_filtered_C3, method='mb', verbose=TRUE, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C3_default) #0.04592922
sum(getRefit(se.mb.breeding_C3_default))/2 #311 #77 #0

betaMat_C3=as.matrix(symBeta(getOptBeta(se.mb.breeding_C3_default)))
positive_C3=length(betaMat_C3[betaMat_C3>0])/2; positive_C3 #34
negative_C3=length(betaMat_C3[betaMat_C3< 0])/2; negative_C3 #0



#se.mb.breeding_C5_0.1_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C5, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C5_0.1_10) #0.04676261
#sum(getRefit(se.mb.breeding_C5_0.1_10))/2 # 364 #59
#se.mb.breeding_C5_0.1_10$select$stars$summary
#se.mb.breeding_C5_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C5, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C5_0.1_30) #0.04583248
#sum(getRefit(se.mb.breeding_C5_0.1_30))/2 # 260 #37
#se.mb.breeding_C5_0.1_30$select$stars$summary
#se.mb.breeding_C5_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C5, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C5_0.01_30) #0.04658554
#sum(getRefit(se.mb.breeding_C5_0.01_30))/2 # 328 #47
#se.mb.breeding_C5_0.01_30$select$stars$summary
#se.mb.breeding_C5_0.01_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C5, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C5_0.01_10) #0.03727863
#sum(getRefit(se.mb.breeding_C5_0.01_10))/2 # 330 #68
#se.mb.breeding_C5_0.01_10$select$stars$summary
se.mb.breeding_C5_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C5, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs)
getStability(se.mb.breeding_C5_0.1_50) #0.04976817
sum(getRefit(se.mb.breeding_C5_0.1_50))/2 # 250 #35
se.mb.breeding_C5_0.1_50$select$stars$summary
se.mb.breeding_C5_default<- spiec.easi(Bacteria_seed_breeding_filtered_C5, method='mb', verbose=TRUE, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C5_default) #0.04456765
sum(getRefit(se.mb.breeding_C5_default))/2 #366 #62 #2

betaMat_C5=as.matrix(symBeta(getOptBeta(se.mb.breeding_C5_default)))
positive_C5=length(betaMat_C5[betaMat_C5>0])/2; positive_C5 #37
negative_C5=length(betaMat_C5[betaMat_C5< 0])/2; negative_C5 #0



#se.mb.breeding_C7_0.1_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C7, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C7_0.1_10) #0.04397511 #0.03003834
#sum(getRefit(se.mb.breeding_C7_0.1_10))/2 #63 #0
#se.mb.breeding_C7_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C7, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C7_0.1_30) #0.04174569
#sum(getRefit(se.mb.breeding_C7_0.1_30))/2 #66
#se.mb.breeding_C7_0.1_30$select$stars$summary #0.01848062           
se.mb.breeding_C7_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C7, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs)
getStability(se.mb.breeding_C7_0.1_50) #0.049559
sum(getRefit(se.mb.breeding_C7_0.1_50))/2 #66
se.mb.breeding_C7_0.1_50$select$stars$summary #0.01848062   
#se.mb.breeding_C7_0.01_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C7, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C7_0.01_50) #0.04537342
#sum(getRefit(se.mb.breeding_C7_0.01_50))/2 #71
#se.mb.breeding_C7_0.01_50$select$stars$summary #0.01912368   
#se.mb.breeding_C7_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C7, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C7_0.01_30) #0.04116312
#sum(getRefit(se.mb.breeding_C7_0.01_30))/2 #74
#se.mb.breeding_C7_0.01_30$select$stars$summary #0.01848062           
#se.mb.breeding_C7_0.02_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C7, method='mb', lambda.min.ratio=0.02, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C7_0.02_30) #0.03821798
#sum(getRefit(se.mb.breeding_C7_0.02_30))/2 #64
#se.mb.breeding_C7_0.02_30$select$stars$summary #0.01848062           
se.mb.breeding_C7_default<- spiec.easi(Bacteria_seed_breeding_filtered_C7, method='mb', verbose=TRUE, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_default) #0.0367634
sum(getRefit(se.mb.breeding_C7_default))/2 #73 #9 #0

betaMat_C7=as.matrix(symBeta(getOptBeta(se.mb.breeding_C7_default)))
positive_C7=length(betaMat_C7[betaMat_C7>0])/2; positive_C7 #0
negative_C7=length(betaMat_C7[betaMat_C7< 0])/2; negative_C7 #0



#se.mb.breeding_C8_0.1_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C8, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C8_0.1_10) #0.04536433
#sum(getRefit(se.mb.breeding_C8_0.1_10))/2 #127 
se.mb.breeding_C8_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C8, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs)
getStability(se.mb.breeding_C8_0.1_30) #0.04746008
sum(getRefit(se.mb.breeding_C8_0.1_30))/2 #100 #0
se.mb.breeding_C8_0.1_30$select$stars$summary #0.01912368   
#se.mb.breeding_C8_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C8, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C8_0.1_50) #0.04684644
#sum(getRefit(se.mb.breeding_C8_0.1_50))/2 #87 #0
#se.mb.breeding_C8_0.1_50$select$stars$summary #0.01912368   
#se.mb.breeding_C8_0.01_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C8, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C8_0.01_10) #0.03697184
#sum(getRefit(se.mb.breeding_C8_0.01_10))/2 #108 #28
#se.mb.breeding_C8_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C8, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs)
#getStability(se.mb.breeding_C8_0.01_30) #0.03657983
se.mb.breeding_C8_0.01_100 <- spiec.easi(Bacteria_seed_breeding_filtered_C8, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=100, pulsar.select='batch', pulsar.params=bargs)
getStability(se.mb.breeding_C8_0.01_100) #0.0467056
sum(getRefit(se.mb.breeding_C8_0.01_100))/2 #187 #0
se.mb.breeding_C8_0.01_100$select$stars$summary #0.01912368   
se.mb.breeding_C8_default<- spiec.easi(Bacteria_seed_breeding_filtered_C8, method='mb', verbose=TRUE, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C8_default) #0.03176934
sum(getRefit(se.mb.breeding_C8_default))/2 #84 #22 #0


betaMat_C8=as.matrix(symBeta(getOptBeta(se.mb.breeding_C8_default)))
positive_C8=length(betaMat_C8[betaMat_C8>0])/2; positive_C8 #24 #0
negative_C8=length(betaMat_C8[betaMat_C8< 0])/2; negative_C8 #0


###BE

#se.mb.breeding_C3_BE_0.1_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C3_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs) #abort
#getStability(se.mb.breeding_C3_BE_0.1_10) #0.03343154
#sum(getRefit(se.mb.breeding_C3_BE_0.1_10))/2 #47 
#se.mb.breeding_C3_BE_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C3_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
#getStability(se.mb.breeding_C3_BE_0.1_30) #0.04311748
#sum(getRefit(se.mb.breeding_C3_BE_0.1_30))/2 #55
#se.mb.breeding_C3_BE_0.1_30$select$stars$summary #0.02018302           
se.mb.breeding_C3_BE_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C3_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C3_BE_0.01_30) #0.04960341
sum(getRefit(se.mb.breeding_C3_BE_0.01_30))/2 #62
se.mb.breeding_C3_BE_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C3_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C3_BE_0.1_50) #0.04278738
sum(getRefit(se.mb.breeding_C3_BE_0.1_50))/2 #53
se.mb.breeding_C3_BE_default<- spiec.easi(Bacteria_seed_breeding_filtered_C3_BE, method='mb', verbose=TRUE, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C3_BE_default) #0.04198943
sum(getRefit(se.mb.breeding_C3_BE_default))/2 #58 #13 #0

betaMat_C3_BE=as.matrix(symBeta(getOptBeta(se.mb.breeding_C3_BE_default)))
positive_C3_BE=length(betaMat_C3_BE[betaMat_C3_BE>0])/2; positive_C3_BE #0
negative_C3_BE=length(betaMat_C3_BE[betaMat_C3_BE< 0])/2; negative_C3_BE #0_BE

getOptMerge(se.mb.breeding_C3_BE_0.01_30) 
sum(getOptMerge(se.mb.breeding_C3_BE_0.01_30) > 0)/2
#https://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php

#se.mb.breeding_C5_BE_0.1_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C5_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs) #abort
#getStability(se.mb.breeding_C5_BE_0.1_10) #0.04135446
#sum(getRefit(se.mb.breeding_C5_BE_0.1_10))/2 #123
#se.mb.breeding_C5_BE_0.1_10$select$stars$summary #0.02018302           
#se.mb.breeding_C5_BE_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C5_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
#getStability(se.mb.breeding_C5_BE_0.1_30) #0.04154749
#sum(getRefit(se.mb.breeding_C5_BE_0.1_30))/2 #102
se.mb.breeding_C5_BE_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C5_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C5_BE_0.01_30) #0.04899283
sum(getRefit(se.mb.breeding_C5_BE_0.01_30))/2 #126
se.mb.breeding_C5_BE_default<- spiec.easi(Bacteria_seed_breeding_filtered_C5_BE, method='mb', verbose=TRUE, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C5_BE_default) #0.04752867
sum(getRefit(se.mb.breeding_C5_BE_default))/2 #141 #34 #0

betaMat_C5_BE=as.matrix(symBeta(getOptBeta(se.mb.breeding_C5_BE_default)))
positive_C5_BE=length(betaMat_C5_BE[betaMat_C5_BE>0])/2; positive_C5_BE #21 (0.1) #23 (0.01)
negative_C5_BE=length(betaMat_C5_BE[betaMat_C5_BE<0])/2; negative_C5_BE #0_BE



se.mb.breeding_C7_BE_0.1_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C7_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_BE_0.1_10) #0.04432807
sum(getRefit(se.mb.breeding_C7_BE_0.1_10))/2 #222 #51
se.mb.breeding_C7_BE_0.1_10$select$stars$summary #0.02018302           
se.mb.breeding_C7_BE_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C7_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_BE_0.1_30) #0.0458316
sum(getRefit(se.mb.breeding_C7_BE_0.1_30))/2 #162 #29
se.mb.breeding_C7_BE_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C7_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_BE_0.01_30) #0.04637208
sum(getRefit(se.mb.breeding_C7_BE_0.01_30))/2 #204 #45
se.mb.breeding_C7_BE_0.01_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C7_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_BE_0.01_50) #0.04748652
sum(getRefit(se.mb.breeding_C7_BE_0.01_50))/2 #175 #34
se.mb.breeding_C7_BE_0.01_50$select$stars$summary #0.02018302           
se.mb.breeding_C7_BE_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C7_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_BE_0.1_50) #0.04706911
sum(getRefit(se.mb.breeding_C7_BE_0.1_50))/2 #132 #0
se.mb.breeding_C7_BE_0.1_100 <- spiec.easi(Bacteria_seed_breeding_filtered_C7_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=100, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_BE_0.1_100) #0.04706911
sum(getRefit(se.mb.breeding_C7_BE_0.1_100))/2 #114 #0
se.mb.breeding_C7_BE_default<- spiec.easi(Bacteria_seed_breeding_filtered_C7_BE, method='mb', verbose=TRUE, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_BE_default) #0.04967014
sum(getRefit(se.mb.breeding_C7_BE_default))/2 #256 #63 #1

betaMat_C7_BE=as.matrix(symBeta(getOptBeta(se.mb.breeding_C7_BE_default)))
positive_C7_BE=length(betaMat_C7_BE[betaMat_C7_BE>0])/2; positive_C7_BE #29
negative_C7_BE=length(betaMat_C7_BE[betaMat_C7_BE< 0])/2; negative_C7_BE #0_BE




se.mb.breeding_C8_BE_0.1_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C8_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C8_BE_0.1_10) #0.04665273
sum(getRefit(se.mb.breeding_C8_BE_0.1_10))/2 #244 #73
se.mb.breeding_C8_BE_0.1_10$select$stars$summary #0.01656383           
se.mb.breeding_C8_BE_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C8_BE, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C8_BE_0.1_30) #0.04491706
sum(getRefit(se.mb.breeding_C8_BE_0.1_30))/2 #186 #29
se.mb.breeding_C8_BE_0.01_10 <- spiec.easi(Bacteria_seed_breeding_filtered_C8_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=10, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C8_BE_0.01_10) #0.04626821
sum(getRefit(se.mb.breeding_C8_BE_0.01_10))/2 #282 #83
se.mb.breeding_C8_BE_0.01_10$select$stars$summary #0.01656383           
se.mb.breeding_C8_BE_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C8_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C8_BE_0.01_30) #0.04285935
sum(getRefit(se.mb.breeding_C8_BE_0.01_30))/2 #282 #47
se.mb.breeding_C8_BE_0.01_10$select$stars$summary #0.01656383           
se.mb.breeding_C8_BE_0.01_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C8_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C8_BE_0.01_50) #0.04893187
sum(getRefit(se.mb.breeding_C8_BE_0.01_50))/2 #214 #44
se.mb.breeding_C8_BE_0.01_50$select$stars$summary #0.01656383           
se.mb.breeding_C8_BE_default<- spiec.easi(Bacteria_seed_breeding_filtered_C8_BE, method='mb', verbose=TRUE, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C8_BE_default) #0.04165089
sum(getRefit(se.mb.breeding_C8_BE_default))/2 #247 #74 #0

betaMat_C8_BE=as.matrix(symBeta(getOptBeta(se.mb.breeding_C8_BE_default)))
positive_C8_BE=length(betaMat_C8_BE[betaMat_C8_BE>0])/2; positive_C8_BE #73
negative_C8_BE=length(betaMat_C8_BE[betaMat_C8_BE<0])/2; negative_C8_BE #0_BE



vsize_all <- log2(apply(n.c_all, 2, mean)) # add log abundance as properties of vertex/nodes.


##final

se.mb.breeding_C3_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C3, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, sel.criterion='stars',pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C3_0.1_50) #0.04794267 
sum(getRefit(se.mb.breeding_C3_0.1_50))/2 # 233 #34

se.mb.breeding_C5_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C5, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs)
getStability(se.mb.breeding_C5_0.1_50) #0.04976817
sum(getRefit(se.mb.breeding_C5_0.1_50))/2 # 250 #35

se.mb.breeding_C7_0.1_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C7, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs)
getStability(se.mb.breeding_C7_0.1_50) #0.049559
sum(getRefit(se.mb.breeding_C7_0.1_50))/2 #66 #0

se.mb.breeding_C8_0.1_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C8, method='mb', lambda.min.ratio=0.1, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs)
getStability(se.mb.breeding_C8_0.1_30) #0.04746008
sum(getRefit(se.mb.breeding_C8_0.1_30))/2 #100 #0

se.mb.breeding_C3_BE_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C3_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C3_BE_0.01_30) #0.04960341
sum(getRefit(se.mb.breeding_C3_BE_0.01_30))/2 #62 #0

se.mb.breeding_C5_BE_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C5_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C5_BE_0.01_30) #0.04899283
sum(getRefit(se.mb.breeding_C5_BE_0.01_30))/2 #126 #23

se.mb.breeding_C7_BE_0.01_30 <- spiec.easi(Bacteria_seed_breeding_filtered_C7_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=30, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C7_BE_0.01_30) #0.04637208
sum(getRefit(se.mb.breeding_C7_BE_0.01_30))/2 #204 #45

se.mb.breeding_C8_BE_0.01_50 <- spiec.easi(Bacteria_seed_breeding_filtered_C8_BE, method='mb', lambda.min.ratio=0.01, verbose=TRUE, nlambda=50, pulsar.select='batch', pulsar.params=bargs) #abort
getStability(se.mb.breeding_C8_BE_0.01_50) #0.04893187
sum(getRefit(se.mb.breeding_C8_BE_0.01_50))/2 #214 #44



ig2.mb_breeding_C3 <- adj2igraph(getRefit(se.mb.breeding_C3_default),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered_C3)))
plot_network(ig2.mb_breeding_C3, Bacteria_seed_breeding_filtered_C3, type='taxa', color = "Phylum", label=NULL)
transitivity(ig2.mb_breeding_C3) #0.2492582

ig2.mb_breeding_C5 <- adj2igraph(getRefit(se.mb.breeding_C5_default),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered_C5)))
plot_network(ig2.mb_breeding_C5, Bacteria_seed_breeding_filtered_C5, type='taxa', color = "Family", label="Genus")
transitivity(ig2.mb_breeding_C5) #0.2982327

ig2.mb_breeding_C7 <- adj2igraph(getRefit(se.mb.breeding_C7_default),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered_C7)))
plot_network(ig2.mb_breeding_C7, Bacteria_seed_breeding_filtered_C7, type='taxa', color = "Family")
transitivity(ig2.mb_breeding_C7) #0

ig2.mb_breeding_C8 <- adj2igraph(getRefit(se.mb.breeding_C8_default),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered_C8)))
plot_network(ig2.mb_breeding_C8, Bacteria_seed_breeding_filtered_C8, type='taxa', color = "Family")
transitivity(ig2.mb_breeding_C8) #0.07594937

ig2.mb_breeding_C3_BE <- adj2igraph(getRefit(se.mb.breeding_C3_BE_default),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered_C3_BE)))
plot_network(ig2.mb_breeding_C3_BE, Bacteria_seed_breeding_filtered_C3_BE, type='taxa', color = "Family")
transitivity(ig2.mb_breeding_C3_BE) #0.2112676

ig2.mb_breeding_C5_BE <- adj2igraph(getRefit(se.mb.breeding_C5_BE_default),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered_C5_BE)))
plot_network(ig2.mb_breeding_C8_BE, Bacteria_seed_breeding_filtered_C8_BE, type='taxa', color = "Family")
transitivity(ig2.mb_breeding_C5_BE) #0.3283019

ig2.mb_breeding_C7_BE <- adj2igraph(getRefit(se.mb.breeding_C7_BE_default),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered_C7_BE)))
plot_network(ig2.mb_breeding_C3_BE, Bacteria_seed_breeding_filtered_C3_BE, type='taxa', color = "Family")
transitivity(ig2.mb_breeding_C7_BE) #0.3024862

ig2.mb_breeding_C8_BE <- adj2igraph(getRefit(se.mb.breeding_C8_BE_default),  vertex.attr=list(name=taxa_names(Bacteria_seed_breeding_filtered_C8_BE)))
plot_network(ig2.mb_breeding_C8_BE, Bacteria_seed_breeding_filtered_C8_BE, type='taxa', color = "Family")
transitivity(ig2.mb_breeding_C8_BE) #0.3352685


betaMat_C3=as.matrix(symBeta(getOptBeta(se.mb.breeding_C8_BE_0.02_500)))
betaMat_C5=as.matrix(symBeta(getOptBeta(se.mb.breeding_C5_0.5_100)))
betaMat_C7=as.matrix(symBeta(getOptBeta(se.mb.breeding_C7_0.5_100)))
betaMat_C8=as.matrix(symBeta(getOptBeta(se.mb.breeding_C8_0.5_100)))

positive_C3=length(betaMat_C3[betaMat_C3>0.3])/2; positive_C3 #0
negative_C3=length(betaMat_C3[betaMat_C3<-0.3])/2; negative_C3 #0


positive_C5=length(betaMat_C5[betaMat_C5>0.3])/2; positive_C5 #86
negative_C5=length(betaMat_C5[betaMat_C5<-0.3])/2; negative_C5 #0

positive_C7=length(betaMat_C7[betaMat_C7>0.3])/2; positive_C7 #0
negative_C7=length(betaMat_C7[betaMat_C7<-0.3])/2; negative_C7 #0

positive_C8=length(betaMat_C8[betaMat_C8>0.3])/2; positive_C8 #24
negative_C8=length(betaMat_C8[betaMat_C8<-0.3])/2; negative_C8 #0

###https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html
install.packages("intergraph")
install.packages("GGally")
devtools::install_github("briatte/ggnet")

install.packages("network")
install.packages("ggnetwork")

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(dplyr) # data handling
library(SpiecEasi) # Network analysis for sparse compositional data  
library(network)
library(intergraph)
#devtools::install_github("briatte/ggnet")
library(ggnet)
library(igraph)
#devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)

network_colors<-c("Acidobacteriota"="#800000", 
                  "Actinobacteriota"="#469990" , 
                  "Bacteroidota"="#808000" , 
                  "Bdellovibrionota"="#9A6324" , 
                  "Campilobacterota"="#000075",
                  "Chloroflexi"="#000000" , 
                  "Cyanobacteria"="#bfef45" , 
                  "Deinococcota"="#4363d8" , 
                  "Desulfobacterota"="#ffe119" ,
                  "Firmicutes"="#e6194B" , 
                  "Fusobacteriota"="#ffd8b1",
                  "Gemmatimonadota"="#42d4f4" ,
                  "Myxococcota"="#3cb44b" , 
                  "NB1-j"="#911eb4" ,
                  "Patescibacteria"="#f032e6",
                  "Planctomycetota"="#fabed4" ,
                  "Proteobacteria"="#f58231" ,
                  "Spirochaetota"="#fffac8",
                  "Synergistota"="#aaffc3",
                  "Verrucomicrobiota"="#dcbeff")  



###C3
Bacteria_seed_breeding_filtered_C3.f <- microbiomeutilities::format_to_besthit(Bacteria_seed_breeding_filtered_C3)
head(tax_table(Bacteria_seed_breeding_filtered_C3.f))
head(tax_table(Bacteria_seed_breeding_filtered_C3.f))
colnames(tax_table(Bacteria_seed_breeding_filtered_C3.f))

OTU3 <- t(otu_table(Bacteria_seed_breeding_filtered_C3.f))


n.c_C3 <- symBeta(getOptBeta(se.mb.breeding_C3_default))
colnames(n.c_C3) <- rownames(n.c_C3) <- colnames(OTU3)


vsize_C3 <- log2(apply(OTU3, 2, mean)) # add log abundance as properties of vertex/nodes.

nw_breeding_C3 <- graph.adjacency(n.c_C3, mode='undirected', add.rownames = TRUE, weighted = TRUE)



family_C3 <- map_levels(colnames(OTU5), from = "best_hit", to = "Family", tax_table(Bacteria_seed_breeding_filtered_C3.f))
phylum_C3 <- map_levels(colnames(OTU5), from = "best_hit", to = "Phylum", tax_table(Bacteria_seed_breeding_filtered_C3.f))
genus_C3 <- map_levels(colnames(OTU5), from = "best_hit", to = "Genus", tax_table(Bacteria_seed_breeding_filtered_C3.f))


E(nw_breeding_C3)$color<-as.character(cut(E(nw_breeding_C3)$weight, 
                                          breaks=c( -1, 0, 1), 
                                          labels=c("red",  "green")))

nw_breeding_C3.net <- asNetwork(nw_breeding_C3)

nw_breeding_C3.net %v% "nodesize" <- vsize_C3
nw_breeding_C3.net %v% "Family" <- family_C3
nw_breeding_C3.net %v% "Phylum" <- phylum_C3
nw_breeding_C3.net %v% "Genus" = genus_C3
#network plot

p3 <- ggnet2(nw_breeding_C3.net, node.color = "Phylum", palette= network_colors,
             label = "Genus", node.size = "nodesize",
             label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p3

p3_unlab <- ggnet2(nw_breeding_C3.net, node.color = "Phylum", palette= network_colors,
                   label = "", node.size = "nodesize",
                   label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p3_unlab



hub_C3<-hub_score(nw_breeding_C3, scale = TRUE, weights = NULL, options = arpack_defaults)
hub_C3<- hub_C3$vector
max(hub_C3)
tail(sort(hub_C3),10)




E(nw_breeding_C5)$color<-as.character(cut(E(nw_breeding_C5)$weight, 
                                          breaks=c( -1, 0, 1), 
                                          labels=c("red",  "green")))


nw_breeding_C5.net <- asNetwork(nw_breeding_C5)

nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)




nw_breeding_C5.net %v% "nodesize" <- vsize_C5
nw_breeding_C5.net %v% "Family" <- family_C5
nw_breeding_C5.net %v% "Phylum" <- phylum_C5
nw_breeding_C5.net %v% "Genus" = genus_C5
#network plot

p5 <- ggnet2(nw_breeding_C5.net, node.color = "Phylum", palette= network_colors,
             label = "Genus", node.size = "nodesize",
             label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p5

p5_unlab <- ggnet2(nw_breeding_C5.net, node.color = "Phylum", palette= network_colors,
                   label = "", node.size = "nodesize",
                   label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p5_unlab


hub_C5<-hub_score(nw_breeding_C5, scale = TRUE, weights = NULL, options = arpack_defaults)
hub_C5<- hub_C5$vector
max(hub_C5)
tail(sort(hub_C5),5)

aut_C5 <- authority_score(nw_breeding_C5, weights=NA)$vector
tail(sort(aut_C5),5)


#####C7
Bacteria_seed_breeding_filtered_C7.f <- microbiomeutilities::format_to_besthit(Bacteria_seed_breeding_filtered_C7)
head(tax_table(Bacteria_seed_breeding_filtered_C7.f))
head(tax_table(Bacteria_seed_breeding_filtered_C7.f))
colnames(tax_table(Bacteria_seed_breeding_filtered_C7.f))

OTU7 <- t(otu_table(Bacteria_seed_breeding_filtered_C7.f))


n.c_C7 <- symBeta(getOptBeta(se.mb.breeding_C7_default))
colnames(n.c_C7) <- rownames(n.c_C7) <- colnames(OTU7)


vsize_C7 <- log2(apply(OTU7, 2, mean)) # add log abundance as properties of vertex/nodes.

nw_breeding_C7 <- graph.adjacency(n.c_C7, mode='undirected', add.rownames = TRUE, weighted = TRUE)



E(nw_breeding_C7)$color<-as.character(cut(E(nw_breeding_C7)$weight, 
                                          breaks=c( -1, 0, 1), 
                                          labels=c("red",  "green")))

nw_breeding_C7.net <- asNetwork(nw_breeding_C7)


family_C7 <- map_levels(colnames(OTU7), from = "best_hit", to = "Family", tax_table(Bacteria_seed_breeding_filtered_C7.f))
phylum_C7 <- map_levels(colnames(OTU7), from = "best_hit", to = "Phylum", tax_table(Bacteria_seed_breeding_filtered_C7.f))
genus_C7 <- map_levels(colnames(OTU7), from = "best_hit", to = "Genus", tax_table(Bacteria_seed_breeding_filtered_C7.f))



nw_breeding_C7.net %v% "nodesize" <- vsize_C7
nw_breeding_C7.net %v% "Family" <- family_C7
nw_breeding_C7.net %v% "Phylum" <- phylum_C7
nw_breeding_C7.net %v% "Genus" = genus_C7
#network plot

p7 <- ggnet2(nw_breeding_C7.net, node.color = "Phylum", palette= network_colors,
             label = "Genus", node.size = "nodesize",
             label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p7

p7_unlab <- ggnet2(nw_breeding_C7.net, node.color = "Phylum", palette= network_colors,
                   label = "", node.size = "nodesize",
                   label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p7_unlab




#####C8
Bacteria_seed_breeding_filtered_C8.f <- microbiomeutilities::format_to_besthit(Bacteria_seed_breeding_filtered_C8)
head(tax_table(Bacteria_seed_breeding_filtered_C8.f))
head(tax_table(Bacteria_seed_breeding_filtered_C8.f))
colnames(tax_table(Bacteria_seed_breeding_filtered_C8.f))

OTU8 <- t(otu_table(Bacteria_seed_breeding_filtered_C8.f))


n.c_C8 <- symBeta(getOptBeta(se.mb.breeding_C8_default))
colnames(n.c_C8) <- rownames(n.c_C8) <- colnames(OTU8)


vsize_C8 <- log2(apply(OTU8, 2, mean)) # add log abundance as properties of vertex/nodes.

nw_breeding_C8 <- graph.adjacency(n.c_C8, mode='undirected', add.rownames = TRUE, weighted = TRUE)



family_C8 <- map_levels(colnames(OTU8), from = "best_hit", to = "Family", tax_table(Bacteria_seed_breeding_filtered_C8.f))
phylum_C8 <- map_levels(colnames(OTU8), from = "best_hit", to = "Phylum", tax_table(Bacteria_seed_breeding_filtered_C8.f))
genus_C8 <- map_levels(colnames(OTU8), from = "best_hit", to = "Genus", tax_table(Bacteria_seed_breeding_filtered_C8.f))


E(nw_breeding_C8)$color<-as.character(cut(E(nw_breeding_C8)$weight, 
                                          breaks=c( -1, 0, 1), 
                                          labels=c("red",  "green")))




nw_breeding_C8.net <- asNetwork(nw_breeding_C8)


nw_breeding_C8.net %v% "nodesize" <- vsize_C8
nw_breeding_C8.net %v% "Family" <- family_C8
nw_breeding_C8.net %v% "Phylum" <- phylum_C8
nw_breeding_C8.net %v% "Genus" = genus_C8
#network plot

p8 <- ggnet2(nw_breeding_C8.net, node.color = "Phylum", palette= network_colors,
             label = "Genus", node.size = "nodesize",
             label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p8


p8_unlab <- ggnet2(nw_breeding_C8.net, node.color = "Phylum", palette= network_colors,
                   label = "", node.size = "nodesize",
                   label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p8_unlab



###BELGIUM
##C3_BE
Bacteria_seed_breeding_filtered_C3_BE.f <- microbiomeutilities::format_to_besthit(Bacteria_seed_breeding_filtered_C3_BE)
head(tax_table(Bacteria_seed_breeding_filtered_C3_BE.f))
head(tax_table(Bacteria_seed_breeding_filtered_C3_BE.f))
colnames(tax_table(Bacteria_seed_breeding_filtered_C3_BE.f))

OTU3_BE <- t(otu_table(Bacteria_seed_breeding_filtered_C3_BE.f))


n.c_C3_BE <- symBeta(getOptBeta(se.mb.breeding_C3_BE_default))
colnames(n.c_C3_BE) <- rownames(n.c_C3_BE) <- colnames(OTU3_BE)


vsize_C3_BE <- log2(apply(OTU3_BE, 2, mean)) # add log abundance as properties of vertex/nodes.

nw_breeding_C3_BE <- graph.adjacency(n.c_C3_BE, mode='undirected', add.rownames = TRUE, weighted = TRUE)



family_C3_BE <- map_levels(colnames(OTU3_BE), from = "best_hit", to = "Family", tax_table(Bacteria_seed_breeding_filtered_C3_BE.f))
phylum_C3_BE <- map_levels(colnames(OTU3_BE), from = "best_hit", to = "Phylum", tax_table(Bacteria_seed_breeding_filtered_C3_BE.f))
genus_C3_BE <- map_levels(colnames(OTU3_BE), from = "best_hit", to = "Genus", tax_table(Bacteria_seed_breeding_filtered_C3_BE.f))


E(nw_breeding_C3_BE)$color<-as.character(cut(E(nw_breeding_C3_BE)$weight, 
                                             breaks=c( -1, 0, 1), 
                                             labels=c("red",  "green")))

nw_breeding_C3_BE.net <- asNetwork(nw_breeding_C3_BE)




nw_breeding_C3_BE.net %v% "nodesize" <- vsize_C3_BE
nw_breeding_C3_BE.net %v% "Family" <- family_C3_BE
nw_breeding_C3_BE.net %v% "Phylum" <- phylum_C3_BE
nw_breeding_C3_BE.net %v% "Genus" = genus_C3_BE
#network plot

p3_BE <- ggnet2(nw_breeding_C3_BE.net, node.color = "Phylum", palette= network_colors,
                label = "Genus", node.size = "nodesize",
                label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p3_BE

p3_BE_unlab <- ggnet2(nw_breeding_C3_BE.net, node.color = "Phylum", palette= network_colors,
                      label = "", node.size = "nodesize",
                      label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p3_BE_unlab



#####C5
Bacteria_seed_breeding_filtered_C5_BE.f <- microbiomeutilities::format_to_besthit(Bacteria_seed_breeding_filtered_C5_BE)
head(tax_table(Bacteria_seed_breeding_filtered_C5_BE.f))
head(tax_table(Bacteria_seed_breeding_filtered_C5_BE.f))
colnames(tax_table(Bacteria_seed_breeding_filtered_C5_BE.f))

OTU5_BE<- t(otu_table(Bacteria_seed_breeding_filtered_C5_BE.f))


n.c_C5_BE <- symBeta(getOptBeta(se.mb.breeding_C5_BE_default))
colnames(n.c_C5_BE) <- rownames(n.c_C5_BE) <- colnames(OTU5_BE)


vsize_C5_BE <- log2(apply(OTU5_BE, 2, mean)) # add log abundance as properties of vertex/nodes.

nw_breeding_C5_BE <- graph.adjacency(n.c_C5_BE, mode='undirected', add.rownames = TRUE, weighted = TRUE)



family_C5_BE <- map_levels(colnames(OTU5_BE), from = "best_hit", to = "Family", tax_table(Bacteria_seed_breeding_filtered_C5_BE.f))
phylum_C5_BE <- map_levels(colnames(OTU5_BE), from = "best_hit", to = "Phylum", tax_table(Bacteria_seed_breeding_filtered_C5_BE.f))
genus_C5_BE <- map_levels(colnames(OTU5_BE), from = "best_hit", to = "Genus", tax_table(Bacteria_seed_breeding_filtered_C5_BE.f))


E(nw_breeding_C5_BE)$color<-as.character(cut(E(nw_breeding_C5_BE)$weight, 
                                             breaks=c( -1, 0, 1), 
                                             labels=c("red",  "green")))

nw_breeding_C5_BE.net <- asNetwork(nw_breeding_C5_BE)


nw_breeding_C5_BE.net %v% "nodesize" <- vsize_C5_BE
nw_breeding_C5_BE.net %v% "Family" <- family_C5_BE
nw_breeding_C5_BE.net %v% "Phylum" <- phylum_C5_BE
nw_breeding_C5_BE.net %v% "Genus" = genus_C5_BE
#network plot

p5_BE <- ggnet2(nw_breeding_C5_BE.net, node.color = "Phylum", palette= network_colors,
                label = "Genus", node.size = "nodesize",
                label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p5_BE
p5_BE_unlab <- ggnet2(nw_breeding_C5_BE.net, node.color = "Phylum", palette= network_colors,
                      label = "", node.size = "nodesize",
                      label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p5_BE_unlab


hub_C5_BE<-hub_score(nw_breeding_C5_BE, scale = TRUE, weights = NULL, options = arpack_defaults)
hub_C5_BE<- hub_C5_BE$vector
max(hub_C5_BE)
tail(sort(hub_C5_BE),5)

aut_C5_BE<-authority_score(nw_breeding_C5_BE, scale = TRUE, weights = NULL, options = arpack_defaults)
aut_C5_BE<- aut_C5_BE$vector
head(sort(aut_C5_BE),5)

#####C7_BE
Bacteria_seed_breeding_filtered_C7_BE.f <- microbiomeutilities::format_to_besthit(Bacteria_seed_breeding_filtered_C7_BE)
head(tax_table(Bacteria_seed_breeding_filtered_C7_BE.f))
head(tax_table(Bacteria_seed_breeding_filtered_C7_BE.f))
colnames(tax_table(Bacteria_seed_breeding_filtered_C7_BE.f))

OTU7_BE <- t(otu_table(Bacteria_seed_breeding_filtered_C7_BE.f))


n.c_C7_BE <- symBeta(getOptBeta(se.mb.breeding_C7_BE_default))
colnames(n.c_C7_BE) <- rownames(n.c_C7_BE) <- colnames(OTU7_BE)


vsize_C7_BE <- log2(apply(OTU7_BE, 2, mean)) # add log abundance as properties of vertex/nodes.

nw_breeding_C7_BE <- graph.adjacency(n.c_C7_BE, mode='undirected', add.rownames = TRUE, weighted = TRUE)



E(nw_breeding_C7_BE)$color<-as.character(cut(E(nw_breeding_C7_BE)$weight, 
                                             breaks=c( -1, 0, 1), 
                                             labels=c("red",  "green")))

nw_breeding_C7_BE.net <- asNetwork(nw_breeding_C7_BE)


family_C7_BE <- map_levels(colnames(OTU7_BE), from = "best_hit", to = "Family", tax_table(Bacteria_seed_breeding_filtered_C7_BE.f))
phylum_C7_BE <- map_levels(colnames(OTU7_BE), from = "best_hit", to = "Phylum", tax_table(Bacteria_seed_breeding_filtered_C7_BE.f))
genus_C7_BE <- map_levels(colnames(OTU7_BE), from = "best_hit", to = "Genus", tax_table(Bacteria_seed_breeding_filtered_C7_BE.f))


nw_breeding_C7_BE.net %v% "nodesize" <- vsize_C7_BE
nw_breeding_C7_BE.net %v% "Family" <- family_C7_BE
nw_breeding_C7_BE.net %v% "Phylum" <- phylum_C7_BE
nw_breeding_C7_BE.net %v% "Genus" = genus_C7_BE
#network plot

p7_BE <- ggnet2(nw_breeding_C7_BE.net, node.color = "Phylum", palette= network_colors,
                label = "Genus", node.size = "nodesize",
                label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p7_BE

p7_BE_unlab <- ggnet2(nw_breeding_C7_BE.net, node.color = "Phylum", palette= network_colors,
                      label = "", node.size = "nodesize",
                      label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p7_BE_unlab


#####C8_BE
Bacteria_seed_breeding_filtered_C8_BE.f <- microbiomeutilities::format_to_besthit(Bacteria_seed_breeding_filtered_C8_BE)
head(tax_table(Bacteria_seed_breeding_filtered_C8_BE.f))
head(tax_table(Bacteria_seed_breeding_filtered_C8_BE.f))
colnames(tax_table(Bacteria_seed_breeding_filtered_C8_BE.f))

OTU8_BE <- t(otu_table(Bacteria_seed_breeding_filtered_C8_BE.f))


n.c_C8_BE <- symBeta(getOptBeta(se.mb.breeding_C8_BE_default))
colnames(n.c_C8_BE) <- rownames(n.c_C8_BE) <- colnames(OTU8_BE)


vsize_C8_BE <- log2(apply(OTU8_BE, 2, mean)) # add log abundance as properties of vertex/nodes.

nw_breeding_C8_BE <- graph.adjacency(n.c_C8_BE, mode='undirected', add.rownames = TRUE, weighted = TRUE)



family_C8_BE <- map_levels(colnames(OTU8_BE), from = "best_hit", to = "Family", tax_table(Bacteria_seed_breeding_filtered_C8_BE.f))
phylum_C8_BE <- map_levels(colnames(OTU8_BE), from = "best_hit", to = "Phylum", tax_table(Bacteria_seed_breeding_filtered_C8_BE.f))
genus_C8_BE <- map_levels(colnames(OTU8_BE), from = "best_hit", to = "Genus", tax_table(Bacteria_seed_breeding_filtered_C8_BE.f))

E(nw_breeding_C8_BE)$color<-as.character(cut(E(nw_breeding_C8_BE)$weight, 
                                             breaks=c( -1, 0, 1), 
                                             labels=c("red",  "green")))


nw_breeding_C8_BE.net <- asNetwork(nw_breeding_C8_BE)


nw_breeding_C8_BE.net %v% "nodesize" <- vsize_C8_BE
nw_breeding_C8_BE.net %v% "Family" <- family_C8_BE
nw_breeding_C8_BE.net %v% "Phylum" <- phylum_C8_BE
nw_breeding_C8_BE.net %v% "Genus" = genus_C8_BE
#network plot

p8_BE <- ggnet2(nw_breeding_C8_BE.net, node.color = "Phylum", palette= network_colors,
                label = "Genus", node.size = "nodesize",
                label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) +ggtitle("Harvested Material C8") ; p8_BE
p8_BE_unlab <- ggnet2(nw_breeding_C8_BE.net, node.color = "Phylum", palette= network_colors,
                      label = "", node.size = "nodesize",
                      label.size = 2, edge.color = "color")+ guides(color=guide_legend(title="Phylum"), size = FALSE) ; p8_BE_unlab



combined_networks<- ggarrange(p3, p5, p7,p8, p3_BE, p5_BE, p7_BE, p8_BE, nrow=2, ncol=4,labels="AUTO", common.legend = TRUE); combined_networks
combined_networks_unlab<- ggarrange(p3_unlab, p5_unlab, p7_unlab,p8_unlab, p3_BE_unlab, p5_BE_unlab, p7_BE_unlab, p8_BE_unlab, nrow=2, ncol=4,labels="AUTO", common.legend = TRUE); combined_networks_unlab


pdf("/mnt/0544f842-c7d0-4f32-b838-5c00df85c63c/Documents/KristinaMichl/NAPERDIV_Samples/combined_networks_unlab.pdf")
