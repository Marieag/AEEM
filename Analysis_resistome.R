#Code written by M. R. Aggerbeck. 
#mrag@envs.au.dk
#Creative commons, 2023. 

# We collected samples once a month from summer to winter.
# The last timepoint was a sample for DNA extraction. 


###----------------PACKAGES--------------
#Packages used in the following code. 
library("phyloseq"); packageVersion("phyloseq")
library("devtools"); packageVersion("devtools")
library("biomformat"); packageVersion("biomformat")
library("vegan"); packageVersion("vegan")
library("ggplot2"); packageVersion("ggplot2")
library("DESeq2"); packageVersion("DESeq2")
library("mvabund"); packageVersion("mvabund")
library("metacoder"); packageVersion("metacoder")
library("taxa"); packageVersion("taxa")
library("microbiomeSeq"); packageVersion("microbiomeSeq")
library("adespatial"); packageVersion("adespatial")
library("ggpubr"); packageVersion("ggpubr")
library("data.table"); packageVersion("data.table")
library("igraph"); packageVersion("igraph")
library("tidyverse"); packageVersion("tidyverse")
library("plotrix"); packageVersion("plotrix")
library("microbiome"); packageVersion("microbiome")
library("ggplotify"); packageVersion("ggplotify")
library("ggordiplots"); packageVersion("ggordiplots")
library("pairwiseAdonis"); packageVersion("pairwiseAdonis")
library("ape"); packageVersion("ape")
library("gridExtra"); packageVersion("gridExtra")
library("directlabels"); packageVersion("directlabels")
###----------------REMEMBER TO RUN FUNCTIONS.R!-----------------

###--------------------Set WD--------------------------

#This page of code is using Phyloseq as the primary package, and uses the others above as additional analyses. 
#Check out the Phyloseq tutorial and the various package manuals for more info on analyses. 

#Set working directory
setwd("O:/Tech_ENVS-EMBI-Afdelingsdrev/Marie/PhD/DATA/RESISTOME_QIIME/Analysis/")
#Set user directory
uzdir <- "O:/Tech_ENVS-EMBI-Afdelingsdrev/Marie/PhD/DATA/RESISTOME_QIIME/Analysis/"

#Specify folder to save current code output in: 

folder = "results_2023_10/"

###----------------Palettes---------------------------------

barplotcol <- c("ivory2", "#d23f56", "#66b747", "#915bc6","#b4b137","#5e79c7", "#db9334","#58acd9","#d6522c","#5bbe7d", "#d2489a","#407f40",
                "#c886ca", "#9b486b", "#4ab29c", "#a45232", "#acaf69", "#db7987", "#957130", "#e1956d", "steelblue1", "gray")

set7 <- c("red", "orange", "gold", "green3", "aquamarine", "royalblue", "darkorchid")
set7dark <- c("red4", "orange4", "gold4", "forestgreen", "aquamarine4", "royalblue4", "darkorchid4")

set6 <- c("tomato", "red", "palegreen", "green3", "steelblue", "royalblue")
set6dark <- c("tomato4", "red4", "palegreen4", "forestgreen", "steelblue4", "royalblue4")

set3 <- c("red", "green3", "royalblue")
set3dark <- c("red4", "forestgreen", "royalblue4")

###----------------Load data--------------------------------


otu_tbl <- read.csv('feature_table.tsv', header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
tax_tbl <- read.csv('tax_table.tsv', header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
meta <- read.csv('metadata.tsv', header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
tree <- read.tree("tree.nwk")

#Clean OTU table
head(otu_tbl)
colnames(otu_tbl) <- gsub("\\.", "_", colnames(otu_tbl))
rownames(otu_tbl) <- otu_tbl$OTU_ID; otu_tbl[1] <- NULL
head(otu_tbl)
dim(otu_tbl)
which(is.na(otu_tbl))

OTU = otu_table(otu_tbl, taxa_are_rows = TRUE)

#Clean TAX table
head(tax_tbl)
colnames(tax_tbl) <- gsub("\\.", "_", colnames(tax_tbl))
rownames(tax_tbl) <- tax_tbl$Feature_ID; tax_tbl[1] <- NULL

tax_tbl <- tax_tbl %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";\\s*")
tax_tbl$Confidence <- NULL
head(tax_tbl)
TAX = tax_table(as.matrix(tax_tbl))

#Clean META table
head(meta)
meta <- as.data.frame(lapply(meta, function(x) gsub("-", "_", x)))

#add columns if necessary
meta <- meta %>%
  mutate(Community = case_when(
    site == "A" ~ "Suburban", site == "D" ~ "Urban", site == "S" ~ "Rural", TRUE ~ NA_character_))

meta <- meta %>%
  mutate(sample_type = case_when(
    site == "A" ~ "sample", site == "D" ~ "sample", site == "S" ~ "sample", 
    site == "Ext" ~ "blank", TRUE ~ NA_character_))


meta <- meta %>%
  mutate(Month = case_when(
    timepoint  == "T1" ~ "July", timepoint  == "T2" ~ "August", timepoint  == "T3" ~ "September",
    timepoint  == "T4" ~ "October", timepoint  == "T5" ~ "November", timepoint  == "T6" ~ "December",
    timepoint  == "T7" ~ "January", TRUE ~ NA_character_))

meta <- meta %>%
  mutate(Sampling_site = case_when(
    site == "A" ~ "Basin", site == "D" ~ "Pipe", site == "S" ~ "Basin", TRUE ~ NA_character_))
#/

rownames(meta) <- meta$SampleID

head(meta)
metadata <- sample_data(meta)

head(OTU); dim(OTU)
head(TAX); dim(TAX)
head(metadata); dim(metadata)
sample_names(metadata)

physeq = phyloseq(OTU, TAX); physeq
physeq <- merge_phyloseq(physeq,metadata,tree)
physeq; sample_names(physeq)

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8194 taxa and 56 samples ]
# sample_data() Sample Data:       [ 56 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 8194 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8194 tips and 8133 internal nodes ]

backup <- physeq

head(otu_table(physeq))
head(tax_table(physeq))
head(sample_data(physeq))

set.seed = 1

###---------------Rarefaction---------------------

physeq_rar = rarefy_even_depth(physeq, 15000, replace=FALSE)
physeq_rar

# 25/7/23
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7566 taxa and 46 samples ]
# sample_data() Sample Data:       [ 46 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 7566 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7566 tips and 7508 internal nodes ]

###---------------Agglomeration------------------

#remember to set tax_glom(NArm) to false
phy_sp <- tax_glom(physeq, taxrank="Species", NArm = FALSE)
phy_sp = prune_taxa(taxa_sums(phy_sp) > 0, phy_sp)
phy_sp 

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1837 taxa and 56 samples ]
# sample_data() Sample Data:       [ 56 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 1837 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1837 tips and 1836 internal nodes ]

phy_sp_rar <- tax_glom(physeq_rar, taxrank="Species", NArm = FALSE)
phy_sp_rar = prune_taxa(taxa_sums(phy_sp_rar) > 0, phy_sp_rar)
phy_sp_rar 

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1767 taxa and 46 samples ]
# sample_data() Sample Data:       [ 46 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 1767 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1767 tips and 1766 internal nodes ]

###-----------Blank removal------------------------------------------

physeq_nb <- subset_samples(physeq, sample_type!="blank")
physeq_nb = prune_taxa(taxa_sums(physeq_nb) > 0, physeq_nb)
physeq_nb
sample_data(physeq_nb)[1:5]
unique(sample_data(physeq_nb)[,11])

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8177 taxa and 54 samples ]
# sample_data() Sample Data:       [ 54 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 8177 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8177 tips and 8116 internal nodes ]

phy_sp_nb <- subset_samples(phy_sp, sample_type!="blank")
phy_sp_nb = prune_taxa(taxa_sums(phy_sp_nb) > 0, phy_sp_nb)
phy_sp_nb
sample_data(phy_sp_nb)[1:5]
unique(sample_data(phy_sp_nb)[,11])

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1835 taxa and 54 samples ]
# sample_data() Sample Data:       [ 54 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 1835 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1835 tips and 1834 internal nodes ]

physeq_rar_nb <- subset_samples(physeq_rar, sample_type!="blank")
physeq_rar_nb = prune_taxa(taxa_sums(physeq_rar_nb) > 0, physeq_rar_nb)
physeq_rar_nb
sample_data(physeq_rar_nb)[1:5]
unique(sample_data(physeq_rar_nb)[,11])

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7566 taxa and 46 samples ]
# sample_data() Sample Data:       [ 46 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 7566 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7566 tips and 7508 internal nodes ]

phy_sp_rar_nb <- subset_samples(phy_sp_rar, sample_type!="blank")
phy_sp_rar_nb = prune_taxa(taxa_sums(phy_sp_rar_nb) > 0, phy_sp_rar_nb)
phy_sp_rar_nb
sample_data(phy_sp_rar_nb)[1:5]
unique(sample_data(phy_sp_rar_nb)[,11])

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1767 taxa and 46 samples ]
# sample_data() Sample Data:       [ 46 samples by 12 sample variables ]
# tax_table()   Taxonomy Table:    [ 1767 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1767 tips and 1766 internal nodes ]

###-----------Make unique names for taxa-----------------------------

physeq
physeq_rar
phy_sp_nb
phy_sp_rar
physeq_rar_nb
phy_sp_rar_nb

temp <- phy_sp_nb
head(tax_table(temp))

mynames = NULL
for (i in 1:length(taxa_names(temp))){
  name <- makeTaxLabel(taxa_names(temp)[i],temp)
  mynames <- rbind(mynames, c(name))
}
#Find duplicates
n_occur <- data.frame(table(mynames)); n_occur[n_occur$Freq > 1,]
taxa_names(temp) = make.unique(mynames, sep="_"); taxa_names(temp)[1:10]
temp; sample_data(temp)[1:5]; otu_table(temp) [1:5]; tax_table(temp)[1:15]

physeq_rar_nb <- temp

###---------------Richness plot-------------------

# Before rarefaction
p <- plot_richness(physeq_nb, x="timepoint", color="site", measures=c("Simpson", "InvSimpson"))
p <- p + geom_boxplot(data = p$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p
head(p$data)
# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="timepoint"
pValueCutoff=0.05
meta_table <- sample_data(physeq_nb)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value and type - usually the second-to-last three.

richness_anova <- p$data[,13:15]
richness_anova <- cbind(richness_anova, p$data$timepoint)
colnames(richness_anova) <- c("samples", "measure", "value", "timepoint")
richness_anova <- richness_anova[order(richness_anova$timepoint),]
head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw #get pairwise p-values
df_pw

#Apply Bonferroni correction
if (is.null(df_pw)) {
  df_pw <- "No significant p-values" 
  file_path <- paste0(folder, "richness_anova", "_output", ".tsv")
  write.table(df_pw, file=file_path, quote=FALSE, sep='\t')
}else {
  if(bonf){  
    numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
    if (is.null(df_pw$p)) {
      bonf.cor <- "No significant p-values"} 
    else {
      anova_res$df_pw <- anova_res$df_pw %>% mutate(
        p = as.factor(as.numeric(as.character(anova_res$df_pw$p)) * numberofsites))
      bonf.cor <- as.numeric(as.matrix(df_pw$p)) * numberofsites}
    temp <- as.factor(bonf.cor)
    df_pw$p <- temp
    df_pw <- df_pw %>% filter(as.numeric(as.character(p)) <= 0.05) %>% mutate(p = as.factor(p))
    file_path <- paste0(folder, "richness_anova", "_bonf_output", ".tsv")
    write.table(df_pw, file=file_path, quote=FALSE, sep='\t')    
  }
}

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)
# p
p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
# p
p<-p+theme_bw()
# p
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
# p
p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")
# p
p<-p+theme(plot.title = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "Alpha Diversity", subtitle = "Overall dataset")

print(p)
print(df_pw)
df_pw <- df_pw %>% 
  filter(as.numeric(as.character(p)) <= 0.05) %>%
  mutate(p = as.factor(p))

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

richness <- p
richness <- richness + scale_color_manual(values = set7dark)+
  scale_fill_manual(values = set7)
print(richness)
file_path <- paste0(folder, "richness", "_output", ".pdf")
ggsave(richness, file=file_path, width = 20, height = 20, units = "cm", useDingbats=F)
print(df_pw)

###===========After rarefaction==================###

#NB! Just for show. Never richness on rarefied data - bias alert! 

p <- plot_richness(physeq_rar_nb, x="timepoint", color="site", measures=c("Simpson", "InvSimpson"))
p <- p + geom_boxplot(data = p$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p

# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="timepoint"
pValueCutoff=0.05
meta_table <- sample_data(phy)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value and type - usually the second-to-last three.

richness_anova <- p$data[,13:15]
richness_anova <- cbind(richness_anova, p$data$timepoint)
colnames(richness_anova) <- c("samples", "measure", "value", "timepoint")
richness_anova <- richness_anova[order(richness_anova$timepoint),]
head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw #get pairwise p-values
df_pw

if (is.null(df_pw)) {
  df_pw <- "No significant p-values" 
  file_path <- paste0(folder, "richness_anova", "_rar_output", ".tsv")
  write.table(df_pw, file=file_path, quote=FALSE, sep='\t')
}else {
  if(bonf){  
    numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
    # bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
    if (is.null(df_pw$p)) {
      bonf.cor <- "No significant p-values"} 
    else {
      bonf.cor <- as.numeric(as.matrix(df_pw$p)) * numberofsites}
    temp <- as.factor(bonf.cor)
    df_pw$p <- temp
    file_path <- paste0(folder, "richness_anova", "_rar_bonf_output", ".tsv")
    write.table(df_pw, file=file_path, quote=FALSE, sep='\t')    
  }
}

#Apply Bonferroni correction


#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)
p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
p<-p+theme_bw()
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")

p<-p+theme(plot.title = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "Symptoms", subtitle = "Overall dataset")

print(p)
print(df_pw)
df_pw <- df_pw %>% 
  filter(as.numeric(as.character(p)) <= 0.05) %>%
  mutate(p = as.factor(p))

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

richness <- p
richness <- richness + scale_color_manual(values = set7dark)+
  scale_fill_manual(values = set7)
print(richness)
file_path <- paste0(folder, "richness", "_rar_output", ".pdf")
ggsave(richness, file=file_path, width = 20, height = 20, units = "cm", useDingbats=F)
print(df_pw)



###---------------Ordination----------------------

# Before rarefaction
phy_sp_nb
pseq.rel <- microbiome::transform(phy_sp_nb, "compositional")
otu <- microbiome::abundances(pseq.rel)
met <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, met$timepoint, type="centroid")
mod

#Calculate eigenvalues and percentage fit
eigenvalues <- mod$eig
total_eigen <- sum(eigenvalues)
percentage_fit <- eigenvalues / total_eigen * 100
percentage_fit

betadisp <- gg_ordiplot(mod, mod$group)
ordination_betadispersion <- betadisp$plot +
  theme(plot.title = element_text(size = 12))+
  labs(title = "Betadispersion", subtitle = "Unrarefied dataset",
       x = paste0("Axis 1 (", round(percentage_fit[1], 2), "%)"),
       y = paste0("Axis 2 (", round(percentage_fit[2], 2), "%)"))+
  theme_bw()+
  theme(legend.title = element_blank())+ 
  scale_color_manual(values = set7)+
  scale_fill_manual(values = set7)+
  theme(strip.background = element_rect(fill="white" )); ordination_betadispersion

# psymp + scale_color_manual(values = set20)+ 
#   scale_fill_manual(values = set20)

file_path <- paste0(folder, "ordination_betadispersion", "_output", ".pdf")
ggsave(ordination_betadispersion, file=file_path, width = 20, height = 20, units = "cm", useDingbats=F)

permanova <- vegan::adonis(t(otu)~timepoint,
                           data = met, permutations=999, method = "bray")
permanova$aov.tab
file_path <- paste0(folder, "permanova", "_output", ".tsv")
write.table(as.data.frame(permanova$aov.tab), file=file_path, quote=FALSE, sep='\t')

post_hoc_permanova <- pairwise.adonis(t(otu), met$timepoint, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova
file_path <- paste0(folder, "post_hoc_permanova", "_output", ".tsv")
write.table(as.data.frame(post_hoc_permanova), file=file_path, quote=FALSE, sep='\t')

###===========After rarefaction==================###

phy_sp_rar
pseq.rel <- microbiome::transform(phy_sp_rar, "compositional")
otu <- microbiome::abundances(pseq.rel)
met <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, met$timepoint, type="centroid")
mod

#Calculate eigenvalues and percentage fit
eigenvalues <- mod$eig
total_eigen <- sum(eigenvalues)
percentage_fit <- eigenvalues / total_eigen * 100
percentage_fit

betadisp <- gg_ordiplot(mod, mod$group)
ordination_betadispersion <- betadisp$plot +
  theme(plot.title = element_text(size = 12))+
  labs(title = "Betadispersion", subtitle = "Rarefied dataset",
       x = paste0("Axis 1 (", round(percentage_fit[1], 2), "%)"),
       y = paste0("Axis 2 (", round(percentage_fit[2], 2), "%)"))+
  theme_bw()+
  theme(legend.title = element_blank())+ 
  scale_color_manual(values = set7)+
  scale_fill_manual(values = set7)+
  theme(strip.background = element_rect(fill="white" )); ordination_betadispersion

file_path <- paste0(folder, "ordination_betadispersion", "_rar_output", ".pdf")
ggsave(ordination_betadispersion, file=file_path, width = 20, height = 20, units = "cm", useDingbats=F)

permanova <- vegan::adonis(t(otu)~timepoint,
                           data = met, permutations=999, method = "bray")
permanova$aov.tab
file_path <- paste0(folder, "permanova", "_rar_output", ".tsv")
write.table(as.data.frame(permanova$aov.tab), file=file_path, quote=FALSE, sep='\t')

post_hoc_permanova <- pairwise.adonis(t(otu), met$timepoint, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova
file_path <- paste0(folder, "post_hoc_permanova", "_rar_output", ".tsv")
write.table(as.data.frame(post_hoc_permanova), file=file_path, quote=FALSE, sep='\t')


###---------------Merge samples for barplots---------------

phy <- phy_sp_rar_nb

sample_data(phy)
samdat <- sample_data(phy)

samdat$site_timepoint <- paste(samdat$site, samdat$timepoint, sep = "_")
samdat
sample_data(phy) <- samdat

unique(sample_data(phy)[,13])

temp = prune_taxa(taxa_sums(phy) > 0, phy)
temp2 <- as.data.frame(unique(sample_data(temp)[,13]))
groups_phy <- as.vector(as.character(temp2$site_timepoint))
sample_data(temp)$group_phy <- get_variable(temp, "site_timepoint") %in% groups_phy

merged_phy = merge_samples(temp, "site_timepoint")

merged_phy
sample_names(temp)
sample_names(merged_phy)

temp <- as.data.frame(sample_data(merged_phy)) 
temp

temp$site <- gsub("1", 'A', temp$site)
temp$site <- gsub("2", 'D', temp$site)
temp$site <- gsub("3", 'S', temp$site)

temp$timepoint <- paste0("T", temp$timepoint)

x <- rownames(temp)
temp$site_timepoint <- x

temp$Community <- ifelse(temp$site == "A", "Suburban",
                                ifelse(temp$site == "D", "Urban",
                                       ifelse(temp$site == "S", "Rural", NA)))

temp$sample_type <- ifelse(temp$site == "A", "sample",
                         ifelse(temp$site == "D", "sample",
                                ifelse(temp$site == "S", "sample", 
                                       ifelse(temp$site == "Ext", "blank", NA))))

temp$Month <- ifelse(temp$timepoint == "T1", "July",
                           ifelse(temp$timepoint == "T2", "August",
                                  ifelse(temp$timepoint == "T3", "September", 
                                         ifelse(temp$timepoint == "T4", "October", 
                                                ifelse(temp$timepoint == "T5", "November", 
                                                       ifelse(temp$timepoint == "T6", "December", 
                                                              ifelse(temp$timepoint == "T7", "January", NA)))))))
temp$Sampling_site <- ifelse(temp$site == "A", "Basin",
                         ifelse(temp$site == "D", "Pipe",
                                ifelse(temp$site == "S", "Basin", NA)))

#/

sample_data(merged_phy) <- temp
sample_data(merged_phy) 

merged_phy
phy

# tmp <- sample_data(merged_phy)
# class(tmp)
# file_path <- paste0(folder, "object", "_output", ".tsv")
# write.table(tmp, file='results_2023/merged_metadata.tsv', quote=FALSE, sep='\t')


###---------------Barplots------------------------

# Before rarefaction

phy_g <- tax_glom(merged_phy, taxrank = "Genus")
head(sample_data(phy_g))

head(otu_table(merged_phy))
# lcbd.rel  = transform_sample_counts(phy_sp, function(x) x / sum(x) )
pseq.rel <- microbiome::transform(merged_phy, "compositional")
head(otu_table(pseq.rel))
p <- plot_taxa_mra(pseq.rel,grouping_column="site_timepoint",method="ab.sorensen")
p
tmp <- as.data.frame(p$data$Sample)
head(tmp)
class(tmp)
tmp <- str_split_fixed(tmp$`p$data$Sample`, "_", 2)
p$data$Site <- tmp[,1]
p$data$Timepoint <- tmp[,2]

head(p$data)
p
p$data$Sample2 <- p$data$Sample
p$data$Sample <- p$data$Groups

# p$data$Sample <- paste0(p$data$Groups, "_", p$data$Sample)
# rownames(p$data) = p$data$Sample
# 
# library("data.table")
# newtab = data.table(p$data)
# setorder(newtab, Groups)
# p$data <- newtab
# print(p)


barplot <- p + 
  labs(title="", subtitle = "")+
  theme(legend.title = element_blank())+
  # scale_color_manual(values = barplotcol)+
  scale_fill_manual(values = barplotcol)#

barplot
head(barplot$data)

barplot_facetgrid <- barplot + facet_wrap(~Timepoint, shrink = TRUE, drop=TRUE, scale="free")+
  labs(y= "Relative Abundance", x = "")

barplot_facetgrid <- barplot_facetgrid + labs(title="20 most common taxa")
barplot_facetgrid
# bp_factgrid2 + theme(strip.text = element_text(size = 12))
file_path <- paste0(folder, "barplot_facetgrid", "_output", ".pdf")
ggsave(barplot_facetgrid, file=file_path, width = 25, height = 15, units = "cm", useDingbats=F)

###===========After rarefaction==================###


# phy_g <- tax_glom(merged_phy, taxrank = "Genus")
head(sample_data(phy_sp_rar))

head(otu_table(phy_sp_rar))
# lcbd.rel  = transform_sample_counts(phy_sp, function(x) x / sum(x) )
pseq.rel <- microbiome::transform(phy_sp_rar, "compositional")
head(otu_table(pseq.rel))
p <- plot_taxa_mra(pseq.rel,grouping_column="timepoint",method="ab.sorensen"); p
tmp <- as.data.frame(p$data$Sample)
head(tmp)
class(tmp)
tmp <- str_split_fixed(tmp$`p$data$Sample`, "_", 2)
p$data$Symptom <- tmp[,1]
tmp <- str_split_fixed(tmp[,2], "_A", 2)
p$data$Cultivar <- tmp[,1]
p$data$Year <- tmp[,2]

head(p$data)
p
p$data$Sample2 <- p$data$Sample
p$data$Sample <- p$data$Groups

# p$data$Sample <- paste0(p$data$Groups, "_", p$data$Sample)
# rownames(p$data) = p$data$Sample
# 
# library("data.table")
# newtab = data.table(p$data)
# setorder(newtab, Groups)
# p$data <- newtab
# print(p)


barplot <- p + 
  labs(title="", subtitle = "")+
  theme(legend.title = element_blank())+
  # scale_color_manual(values = barplotcol)+
  scale_fill_manual(values = barplotcol)#

barplot
head(barplot$data)

barplot_facetgrid <- barplot + facet_grid(~Cultivar+Year, shrink = TRUE, drop=TRUE, scale="free")+
  labs(y= "Relative Abundance", x = "")

barplot_facetgrid <- barplot_facetgrid + labs(title="20 most common taxa")
# bp_factgrid2 + theme(strip.text = element_text(size = 12))
file_path <- paste0(folder, "barplot_facetgrid", "_output", ".pdf")
ggsave(barplot_facetgrid, file=file_path, width = 25, height = 15, units = "cm", useDingbats=F)

###--------------------Analysis of T7-----------------------------------

# Given the difference in DNA quality of the various samples, and that only T4 is of a quality comparable with T7
# we've decided to perform community analysis on these timepoint only. 
# This means the analysis from here on focuses on T4+T7 only. 

# T7 <- subset_samples(phy_sp_nb, timepoint == "T7")
# T7 = prune_taxa(taxa_sums(T7) > 0, T7)
# T7
# 
# 
# 
# T7_rar <- subset_samples(phy_sp_rar_nb, timepoint == "T7")
# T7_rar = prune_taxa(taxa_sums(T7_rar) > 0, T7_rar)
# T7_rar

T7  <- subset_samples(phy_sp_nb, timepoint %in% c("T7"))
T7  = prune_taxa(taxa_sums(T7 ) > 0, T7 )
T7
head(tax_table((T7)))
tax_table(T7)[1:15]
sample_data(T7)

tax_table(T7) <-  gsub("\\bs__uncultured\\b", "", tax_table(T7));
# tax_table(T7) <-  gsub("^_.*", "", tax_table(T7));
tax_table(T7) <-  gsub("s__metagenome", "", tax_table(T7));
tax_table(T7) <-  gsub("s__uncultured_bacterium", "", tax_table(T7));
tax_table(T7) <-  gsub("g__uncultured", "", tax_table(T7)); 
tax_table(T7) <-  gsub("f__uncultured", "", tax_table(T7));
tax_table(T7) <-  gsub("o__uncultured", "", tax_table(T7));
tax_table(T7) <-  gsub("c__uncultured", "", tax_table(T7));
tax_table(T7) <-  gsub("p__uncultured", "", tax_table(T7))
head(tax_table((T7)))


temp <- T7
head(tax_table(temp))

mynames = NULL
for (i in 1:length(taxa_names(temp))){
  name <- makeTaxLabel(taxa_names(temp)[i],temp)
  mynames <- rbind(mynames, c(name))
}
#Find duplicates
n_occur <- data.frame(table(mynames)); n_occur[n_occur$Freq > 1,]
taxa_names(temp) = make.unique(mynames, sep="_"); taxa_names(temp)[1:10]
temp; sample_data(temp)[1:5]; otu_table(temp) [1:5]; tax_table(temp)[1:15]

T7 <- temp

temp <- data.frame(sample_data(T7))
temp <- temp %>%
  mutate(Season = ifelse(Month == "October", "Autumn", "Winter"))
temp <- temp %>%
  mutate(Avg_Temp_C = ifelse(Month == "October", "11.2", "0.7"))
sample_data(T7) <- temp
sample_data(T7)


T7_rar <- subset_samples(phy_sp_rar_nb, timepoint == "T7")
T7_rar = prune_taxa(taxa_sums(T7_rar) > 0, T7_rar)
T7_rar

###---------------Agglomerate and rename-------------------



# tax_table(T7) <-  gsub("s__uncultured", "", tax_table(T7));
# tax_table(T7) <-  gsub("^_.*", "", tax_table(T7));
# tax_table(T7) <-  gsub("s__metagenome", "", tax_table(T7));
# tax_table(T7) <-  gsub("s__uncultured_bacterium", "", tax_table(T7));
# tax_table(T7) <-  gsub("g__uncultured", "", tax_table(T7)); 
# tax_table(T7) <-  gsub("f__uncultured", "", tax_table(T7));
# tax_table(T7) <-  gsub("o__uncultured", "", tax_table(T7));
# tax_table(T7) <-  gsub("c__uncultured", "", tax_table(T7));
# tax_table(T7) <-  gsub("p__uncultured", "", tax_table(T7))
# head(tax_table((T7)))

T7_g <- tax_glom(T7, taxrank = "Genus", NArm = FALSE)

temp <- T7_g
head(tax_table(temp))

mynames = NULL
for (i in 1:length(taxa_names(temp))){
  name <- makeTaxLabel(taxa_names(temp)[i],temp)
  mynames <- rbind(mynames, c(name))
}
#Find duplicates
n_occur <- data.frame(table(mynames)); n_occur[n_occur$Freq > 1,]
taxa_names(temp) = make.unique(mynames, sep="_"); taxa_names(temp)[1:10]
temp; sample_data(temp)[1:5]; otu_table(temp) [1:5]; tax_table(temp)[1:15]

T7_g <- temp

# temp <- data.frame(sample_data(T7_g))
# temp$site_season <- paste(temp$site, temp$Season, sep = "_")
# sample_data(T7_g) <- temp
# sample_data(T7_g)

###---------------Richness plot-------------------

# T7 only
p <- plot_richness(T7_g, x="site", color="site", measures=c("Simpson", "InvSimpson"))
p <- p + geom_boxplot(data = p$data, aes(color = NULL), alpha = 0.1)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(legend.position = "none")+
  theme(strip.background = element_rect(fill="white" )); p
head(p$data)
# method = c("richness","simpson", "evenness")
bonf=TRUE
grouping_column ="site"
pValueCutoff=0.05
meta_table <- sample_data(T7_g)

#create a starting matrix from a previous run - check proper columns are chosen 
#get sample id, variable, value - usually the second-to-last three.

richness_anova <- p$data[,15:17]
richness_anova <- cbind(richness_anova, p$data$site)
colnames(richness_anova) <- c("samples", "measure", "value", "site")
# richness_anova <- richness_anova[order(richness_anova$site),]
head(richness_anova)

#perform anova of diversity measure between groups
anova_res <- perform_anova_mra(richness_anova,meta_table,grouping_column,pValueCutoff)
df_pw <- anova_res$df_pw #get pairwise p-values
df_pw

#Apply Bonferroni correction
if (is.null(df_pw)) {
  df_pw <- "No significant p-values" 
  file_path <- paste0(folder, "richness_anova", "_T7_g", ".tsv")
  write.table(df_pw, file=file_path, quote=FALSE, sep='\t')
}else {
  if(bonf){  
    numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
    if (is.null(df_pw$p)) {
      bonf.cor <- "No significant p-values"} 
    else {
      anova_res$df_pw <- anova_res$df_pw %>% mutate(
        p = as.factor(as.numeric(as.character(anova_res$df_pw$p)) * numberofsites))
      bonf.cor <- as.numeric(as.matrix(df_pw$p)) * numberofsites}
    temp <- as.factor(bonf.cor)
    df_pw$p <- temp
    df_pw <- df_pw %>% filter(as.numeric(as.character(p)) <= 0.05) %>% mutate(p = as.factor(p))
    file_path <- paste0(folder, "richness_anova", "_bonf_T7_g", ".tsv")
    write.table(df_pw, file=file_path, quote=FALSE, sep='\t')    
  }
}

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column, fill=grouping_column),data=richness_anova)
p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
p<-p+theme_bw()
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")
p<-p+theme(plot.title = element_text(size = 12)) +
  theme(strip.background = element_rect(fill = "white"))+xlab("")+ 
  labs(title = "Alpha Diversity", subtitle = "Overall dataset")

print(p)
print(df_pw)
df_pw <- df_pw %>% 
  filter(as.numeric(as.character(p)) <= 0.05) %>%
  mutate(p = as.factor(p))

#This loop will generate the lines and significances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(
    x = c(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"])),
          which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"]))), 
    y = c(as.numeric(as.character(df_pw[i,"y"])),
          as.numeric(as.character(df_pw[i,"y"]))), 
    measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
    color="black",
    lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  p
  p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(
    x=(which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"from"]))+
         which(unique(richness_anova[,grouping_column])==as.character(df_pw[i,"to"])))/2,
    y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
    label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                           breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}

richness <- p
richness <- richness + scale_color_manual(values = set3dark)+
  scale_fill_manual(values = set3)
print(richness)
file_path <- paste0(folder, "richness", "_T7_g", ".pdf")
ggsave(richness, file=file_path, width = 20, height = 20, units = "cm", useDingbats=F)
print(df_pw)

###---------------Ordination----------------------

# Before rarefaction
T7_g
pseq.rel <- microbiome::transform(T7_g, "compositional")
otu <- microbiome::abundances(pseq.rel)
met <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, met$site, type="centroid")
mod

#Calculate eigenvalues and percentage fit
eigenvalues <- mod$eig
total_eigen <- sum(eigenvalues)
percentage_fit <- eigenvalues / total_eigen * 100
percentage_fit
mod$group


betadisp <- gg_ordiplot(mod, mod$group)
ordination_betadispersion <- betadisp$plot +
  theme(plot.title = element_text(size = 12))+
  labs(title = "Betadispersion", subtitle = "",
       x = paste0("Axis 1 (", round(percentage_fit[1], 2), "%)"),
       y = paste0("Axis 2 (", round(percentage_fit[2], 2), "%)"))+
  theme_bw()+
  theme(legend.title = element_blank())+ 
  scale_color_manual(values = set3)+
  scale_fill_manual(values = set3)+
  theme(strip.background = element_rect(fill="white" )); ordination_betadispersion


file_path <- paste0(folder, "ordination_betadispersion", "_T7_g", ".pdf")
ggsave(ordination_betadispersion, file=file_path, width = 20, height = 20, units = "cm", useDingbats=F)

permanova <- vegan::adonis(t(otu)~site,
                           data = met, permutations=999, method = "bray")
permanova$aov.tab
file_path <- paste0(folder, "permanova", "_T7_g", ".tsv")
write.table(as.data.frame(permanova$aov.tab), file=file_path, quote=FALSE, sep='\t')

post_hoc_permanova <- pairwise.adonis(t(otu), met$site, sim.function = "vegdist",
                                      sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                      perm = 999)
post_hoc_permanova
file_path <- paste0(folder, "post_hoc_permanova", "_T7_g", ".tsv")
write.table(as.data.frame(post_hoc_permanova), file=file_path, quote=FALSE, sep='\t')

###---------------Barplot, T7------------------------

head(sample_data(T7_g))
head(otu_table(T7_g))
pseq.rel <- microbiome::transform(T7_g, "compositional")
head(otu_table(pseq.rel))
p <- plot_taxa_mra(pseq.rel,grouping_column="site",method="ab.sorensen")
p
tmp <- as.data.frame(p$data$Sample)
head(tmp)
class(tmp)
tmp <- str_split_fixed(tmp$`p$data$Sample`, "_", 2)
p$data$Site <- tmp[,1]
p$data$Timepoint <- tmp[,2]

head(p$data)
p
# p$data$Sample2 <- p$data$Sample
# p$data$Sample <- p$data$Groups

barplot <- p + 
  labs(title="", subtitle = "")+
  theme(legend.title = element_blank())+
  # scale_color_manual(values = barplotcol)+
  scale_fill_manual(values = barplotcol)#

barplot
head(barplot$data)

barplot_facetgrid <- barplot + facet_wrap(~Site, shrink = TRUE, drop=TRUE, scale="free")+
  labs(y= "Relative Abundance", x = "")

barplot_facetgrid <- barplot_facetgrid + labs(title="20 most common taxa")
barplot_facetgrid
# bp_factgrid2 + theme(strip.text = element_text(size = 12))
file_path <- paste0(folder, "barplot_facetgrid", "_output", ".pdf")
ggsave(barplot_facetgrid, file=file_path, width = 25, height = 15, units = "cm", useDingbats=F)

###------------Heat tree----------------

minTotRelAbun = 1 
x = taxa_sums(T7_g)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
T7_g_1 = prune_taxa(keepTaxa, T7_g)
T7_g_1

obj <- parse_phyloseq(T7_g_1)
tissuegroup <- obj$data$sample_data$site
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id, groups = tissuegroup)
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
hist(obj$data$diff_table$wilcox_p_value)

obj$data$diff_table$wilcox_p_value

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of taxa",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'results_2023_10/heattree_T7_g_1.pdf')

###-----------Test network-------------------

#Filters out anything below 0.1% of total abundance.
minTotRelAbun = 0.1 
x = taxa_sums(T7_g)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
T7_g_01 = prune_taxa(keepTaxa, T7_g)
T7_g_01

unique(tax_table(T7_g_01)[,2])

T7_g_01 <- prune_taxa(!taxa_names(T7_g_01) %in% c(
  'd__Bacteria'), T7_g_01)

taxa_names(T7_g_01)
tax_table(T7_g_01)
otu_table(T7_g_01)


plot_net(T7_g_01, maxdist = 0.3, color = "site", shape="site")
head(tax_table(T7_g))
head(otu_table(T7_g))
head(sample_data(T7_g))

T7_f_01 <- tax_glom(T7_g_01, taxrank = "Family", NArm = FALSE)
T7_f_01

tax_table(T7_f_01)[is.na(tax_table(T7_f_01))] <- ""

temp <- T7_f_01

mynames = NULL
for (i in 1:length(taxa_names(temp))){
  name <- makeTaxLabel(taxa_names(temp)[i],temp)
  mynames <- rbind(mynames, c(name))
}
#Find duplicates
n_occur <- data.frame(table(mynames)); n_occur[n_occur$Freq > 1,]
taxa_names(temp) = make.unique(mynames, sep="_"); taxa_names(temp)[1:10]
temp; sample_data(temp)[1:5]; otu_table(temp) [1:5]; tax_table(temp)[1:15]

T7_f_01 <- temp

set.seed(1)

ig <- make_network(T7_f_01, type="samples", max.dist=0.4)
ig_plot <- plot_network(ig, physeq=T7_f_01, type="samples", 
             color="site", shape="timepoint", point_size=4, alpha=1,
             label="value", hjust = 1.35, 
             line_weight=0.5, line_alpha=0.4,
             layout.method=layout.fruchterman.reingold, title=NULL)
ig_plot

ig2 <- make_network(T7_f_01, type="taxa", max.dist=0.4)
set.seed(1)
ig_plot1 <- plot_network(ig2, physeq=T7_f_01, type="taxa", 
             color="Phylum", point_size=4, alpha=1,
             label="value", hjust = 1.35, 
             line_weight=0.5, line_alpha=0.4,
             layout.method=layout.fruchterman.reingold, title=NULL)
ig_plot1

###-----------Ordisurf, metabolite import----------------------

#Import metabolites
# otu_metabolites <- read.csv('otu_metabolites.tsv', header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
load("otu_metabolites_t6.RData")
head(otu_metabolites)
T7_g
rownames(otu_metabolites)
sample_names(T7_g)

sample_names_ordi <- c("A_s1", "A_s2", "A_s3", "D_s1", "D_s2", "D_s3", "S_s1", "S_s2", "S_s3")

temp <- otu_metabolites
rownames(temp) = sample_names_ordi
head(temp)

ps_ordi <- T7_g_1
sample_names(ps_ordi) = sample_names_ordi
sample_names(ps_ordi)

temp2 <- sample_data(cbind(as.data.frame(sample_data(ps_ordi)), temp))
temp2

ps_ordi@sam_data <- temp2
head(sample_data(ps_ordi))
###-----------Ordisurf--------------------------

#----Pick an object to run ordisurf+plots on---- 

result <- generate_filenames(ps_ordi, "ordi"); ordi

# Access the results
new_obj_name <- result$new_obj_name; new_obj_name
original_name <- result$original_name; original_name
filepath <- result$pvalues_file; filepath
filepath2 <- result$ordisurf_file; filepath2
filepath3 <- paste0(original_name, "_ordisurf_contour.pdf"); filepath3

# ordi <- prune_taxa(!taxa_names(ordi) %in% c(
#   'k__Fungi', 'p__Basidiomycota', 'p__Ascomycota', 'f__Pleosporaceae'), ordi)
# ordi = prune_taxa(taxa_sums(ordi) > 0, ordi); ordi

sample_data(ordi)
tax_table(ordi)

#Correct filenames, or they will be overwritten!!!
# ordi <- physeq_antho_21_RA01
# filepath <- "Anthocyanins_RA01_pvalues.txt"
# filepath2 <- "Anthocyanin_ordisurf_RA01.pdf"
cat(original_name, "\n\n", file = filepath, append = FALSE)

pseq.rel <- microbiome::transform(ordi, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

dist <- vegan::vegdist(t(otu), method="bray")
mod <- vegan::betadisper(dist, meta$site, type="centroid")
mod

#Calculate eigenvalues and percentage fit
eigenvalues <- mod$eig
total_eigen <- sum(eigenvalues)
percentage_fit <- eigenvalues / total_eigen * 100
percentage_fit

#Update axis labels
axis_labels <- paste0("Axis ", 1:length(eigenvalues), " (", round(percentage_fit, 2), "%)")


###----------Ordisurf-------------

##https://chrischizinski.github.io/rstats/ordisurf/ - check this to convert to ggplot!!! 

#Ordisurf otu table

set.seed(1)


identical(rownames(meta),rownames(t(otu)))
# temp <- make.unique(as.character(meta$Symptom_Year_Cultivar), sep="_")
# temp
# rownames(meta) <- temp
mydata <- t(otu)
# rownames(mydata) <- temp
meta; dim(meta)
mydata

dist.f <- vegdist(mydata, method = "bray")
# dist.f2  <- stepacross(dist.f, path = "extended")

NMDS <- metaMDS(mydata, noshare = FALSE, autotransform = FALSE, trymax = 500)
# NMDS <- metaMDS(dist.f2, trymax = 500)
sppscores(NMDS) <- mydata
NMDS

plot(NMDS, type="n")
points(NMDS, display = "sites", cex = 1, pch = 16, col = "red")
points(NMDS, display = "species", cex = 1, col = "blue")
# 

# results<-metaMDS(dist.f2, trymax=500)

###----GGPLOT ORDISURF-------------------

# filepath <- "ordisurf_pvalues.txt"

#Total
###ggplot ordisurf
species.scores <- as.data.frame(scores(NMDS, "species"))
species.scores$species <- rownames(species.scores)
names(species.scores)[c(1, 2)] <- c("x", "y")
species.scores$z <- NA
head(species.scores)

sites.scores <- as.data.frame(scores(NMDS, "sites"))
sites.scores$sites <- rownames(sites.scores)
sites.scores <- separate(sites.scores, sites, into = c("Site", "Sample"), sep = "_")
sites.scores$sites <- rownames(sites.scores)
names(sites.scores)[c(1, 2)] <- c("x", "y")
sites.scores$z <- NA
head(sites.scores)

# head(meta)

dist.sf <- ordisurf(NMDS ~ X_Erythromycin  , data = meta, knots= 3, plot = FALSE)
summary(dist.sf)
cat("Personal Care Products", "\n", file = filepath, append = TRUE)
summary_output <- capture.output(summary(dist.sf))
summary_output <- paste(summary_output, collapse = "\n")
cat(summary_output, "\n\n", file = filepath, append = TRUE)
p1 <- plot(dist.sf, col = "forestgreen", add = TRUE)
p1

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}

contour_ordi <- extract.xyz(obj = dist.sf)
head(contour_ordi)
p_start <- ggplot(data = contour_ordi, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1.5)) 
p_start <- direct.label(p_start, "bottom.pieces")
print(p_start)
# class(sites.scores$Symptom)
# class(sites.scores$Year)

# Simplified labels
label_map_test <- data.frame(
  species = unique(species.scores$species),
  label_species = sprintf("%02d", 1:length(unique(species.scores$species)))
); label_map_test

species_scores_test <- species.scores %>%
  left_join(label_map_test, by = "species"); species_scores_test


p_key <- ggplot()  + 
  geom_point(data = sites.scores, aes(x = x, y = y, shape = Site, color = Site), size = 3) +
  # geom_point(data = species.scores, aes(x = x, y = y), shape = 16, size = 1, color = "blue" ) +
  geom_text(data = species_scores_test, check_overlap = FALSE, colour = "blue", aes(
    x = x, y = y, label = label_species), size = 3) +
  theme_bw()+
  coord_equal() + 
  labs(title = "Taxa", 
       x = paste0("NMDS 1 (", round(percentage_fit[1], 2), "%)"),
       y = paste0("NMDS 2 (", round(percentage_fit[2], 2), "%)")) +
  # scale_colour_manual(values = c("darkorchid1", "red3")) +
  theme(legend.position = "bottom")#, legend.box = "vertical")
print(p_key)

# p_wrap <- p_key +
#   facet_wrap(~ Cultivar) +
#   theme(strip.background = element_rect(fill = "white", color = NA),
#         strip.text = element_text(color = "black", face = "bold"))+
#   theme(legend.position = "none"); p_wrap

# p_key_both <- ggplot()  + 
#   geom_point(data = sites.scores, aes(x = x, y = y, shape = Cultivar, color = Cultivar), size = 3) +
#   geom_point(data = species_scores_test, aes(x = x, y = y), shape = 16, size = 1, color = "blue" ) +
#   geom_text(data = species_scores_test, check_overlap = FALSE, colour = "blue", aes(
#     x = x, y = y, label = label_species), size = 3); print(p_key_both)

#legend 
# Create a lookup table for the legend
legend_table_test <- label_map_test %>%
  arrange(label_species); legend_table_test

# Remove row names
rownames(legend_table_test) <- legend_table_test$label_species

legend_table_test <- legend_table_test %>%
  arrange(label_species) %>%
  select(-label_species); legend_table_test

legend_table_test <- legend_table_test %>%
  rename(Species = species) %>%
  mutate(Species = str_replace_all(Species, "o__", "(o) "),
         Species = str_replace_all(Species, "f__", "(f) "),
         Species = str_replace_all(Species, "g__", "(g) "),
         Species = str_replace_all(Species, "s__", "(s) "),
         Species = str_replace_all(Species, "_", " "))

print(legend_table_test)

# Convert the legend table to a graphical object and style it
theme_custom <- ttheme_default(
  core = list(bg_params = list(fill = "white", col = "black"),
              fg_params = list(fontface = 1, cex = 0.8)),
  colhead = list(bg_params = list(fill = "white", col = NA),
                 fg_params = list(fontface = 2, cex = 0.8)),
  rowhead = list(bg_params = list(fill = "white", col = "white"),
                 fg_params = list(fontface = 2, cex = 0.8, col = "blue"))
)


legend_grob_test <- tableGrob(legend_table_test, theme = theme_custom); legend_grob_test

# Arrange the plot and the legend side by side
grid.arrange(legend_grob_test, p_key, ncol = 2, widths = c(1,2))


p_ordi <- p_start + 
  geom_point(data = sites.scores, aes(x = x, y = y, shape = Site), size = 3 ) +
  # geom_point(data = species.scores, aes(x = x, y = y), shape = 16, size = 1, color = "blue" ) +
  geom_text(data = species_scores_test, check_overlap = FALSE, colour = "blue", aes(
    x = x, y = y, label = label_species), size = 3) +
  theme_bw()+
  theme(legend.position = "bottom")+
  coord_equal() +
  labs(title = "Total Industrial Chemicals", 
       x = paste0("NMDS 1 (", round(percentage_fit[1], 2), "%)"),
       y = paste0("NMDS 2 (", round(percentage_fit[2], 2), "%)")) 
print(p_ordi)
p_ordi_indchem <- p_ordi

p_ordi_pcpp
p_ordi_indchem
p_ordi_ther
p_ordi_anti
p_ordi_erythr
p_ordi_pesti

###-----Tri------------

plot(NMDS, type="n")
points(NMDS, display = "sites", cex = 1, pch = 16, col = "red")
points(NMDS, display = "species", cex = 1, col = "blue")
dist.sf <- ordisurf(NMDS ~ Tri_hydroxylated_imputed, data = meta, plot = FALSE)
summary(dist.sf)
cat("Tri hydroxylated", "\n", file = filepath, append = TRUE)
summary_output <- capture.output(summary(dist.sf))
summary_output <- paste(summary_output, collapse = "\n")
cat(summary_output,"\n\n", file = filepath, append = TRUE)
p1 <- plot(dist.sf, col = "forestgreen", add = TRUE)
p1

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}

contour_tri <- extract.xyz(obj = dist.sf)
head(contour_tri)
p_start <- ggplot(data = contour_tri, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1.5)) 
p_start <- direct.label(p_start, "bottom.pieces")
print(p_start)

p_ordi <- p_start + 
  geom_point(data = sites.scores, aes(x = x, y = y, shape = Cultivar), size = 3 ) +
  geom_point(data = species.scores, aes(x = x, y = y), shape = 16, size = 1, color = "blue" ) +
  theme_bw()+
  theme(legend.position = "none")+
  coord_equal() + 
  labs(x = paste0("NMDS 1 (", round(percentage_fit[1], 2), "%)"),
       y = paste0("NMDS 2 (", round(percentage_fit[2], 2), "%)"), 
       title = "Tri-hydroxylated")
print(p_ordi)

p_ordi_tri <- p_ordi +
  facet_wrap(~ Cultivar) +
  theme(strip.background = element_rect(fill = "white", color = NA),
        strip.text = element_text(color = "black", face = "bold"))


###-----Di-------------

plot(NMDS, type="n")
points(NMDS, display = "sites", cex = 1, pch = 16, col = "red")
points(NMDS, display = "species", cex = 1, col = "blue")
dist.sf <- ordisurf(NMDS ~ Di_hydroxylated_imputed, data = meta, plot = FALSE)
summary(dist.sf)
cat("Di hydroxylated", "\n", file = filepath, append = TRUE)
summary_output <- capture.output(summary(dist.sf))
summary_output <- paste(summary_output, collapse = "\n")
cat(summary_output,"\n\n", file = filepath, append = TRUE)
p1 <- plot(dist.sf, col = "forestgreen", add = TRUE)
p1

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}

contour_di <- extract.xyz(obj = dist.sf)
head(contour_di)
p_start <- ggplot(data = contour_di, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1.5)) 
p_start <- direct.label(p_start, "bottom.pieces")
print(p_start)

p_ordi <- p_start + 
  geom_point(data = sites.scores, aes(x = x, y = y, shape = Cultivar), size = 3 ) +
  geom_point(data = species.scores, aes(x = x, y = y), shape = 16, size = 1, color = "blue" ) +
  theme_bw()+
  theme(legend.position = "none")+
  coord_equal() + 
  labs(x = paste0("NMDS 1 (", round(percentage_fit[1], 2), "%)"),
       y = paste0("NMDS 2 (", round(percentage_fit[2], 2), "%)"), 
       title = "Di-hydroxylated")
print(p_ordi)


p_ordi_di <- p_ordi+
  facet_wrap(~ Cultivar) +
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black", face = "bold"))


###-----Acyl-----------

plot(NMDS, type="n")
points(NMDS, display = "sites", cex = 1, pch = 16, col = "red")
points(NMDS, display = "species", cex = 1, col = "blue")
dist.sf <- ordisurf(NMDS ~ Acyl_derivatives_imputed, data = meta, plot = FALSE)
summary(dist.sf)
cat("Acyl derivatives", "\n", file = filepath, append = TRUE)
summary_output <- capture.output(summary(dist.sf))
summary_output <- paste(summary_output, collapse = "\n")
cat(summary_output,"\n\n", file = filepath, append = TRUE)

p1 <- plot(dist.sf, col = "forestgreen", add = TRUE)
p1

extract.xyz <- function(obj) {
  xy <- expand.grid(x = obj$grid$x, y = obj$grid$y)
  xyz <- cbind(xy, c(obj$grid$z))
  names(xyz) <- c("x", "y", "z")
  return(xyz)
}

contour_acyl <- extract.xyz(obj = dist.sf)
head(contour_acyl)
p_start <- ggplot(data = contour_acyl, aes(x, y, z = z)) + stat_contour(aes(colour = ..level..)) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1.5)) 
p_start <- direct.label(p_start, "bottom.pieces")
print(p_start)

p_ordi <- p_start + 
  geom_point(data = sites.scores, aes(x = x, y = y, shape = Cultivar), size = 3 ) +
  geom_point(data = species.scores, aes(x = x, y = y), shape = 16, size = 1, color = "blue" ) +
  theme_bw()+
  theme(legend.position = "none")+
  coord_equal() + 
  labs(x = paste0("NMDS 1 (", round(percentage_fit[1], 2), "%)"),
       y = paste0("NMDS 2 (", round(percentage_fit[2], 2), "%)"), 
       title = "Acyl derivatives")
print(p_ordi)


p_ordi_acyl <- p_ordi +
  facet_wrap(~ Cultivar) +
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(color = "black", face = "bold"))

###----Combine anthocyanin plots-----------


p_key
p_ordi_pcpp
p_ordi_indchem
p_ordi_ther
p_ordi_anti
p_ordi_erythr
p_ordi_pesti

anti_combo <- ggarrange(p_ordi_anti, p_ordi_erythr, ncol = 1, nrow = 2)
anti_combo

# P_wrap3 <- ggarrange(p_wrap, p_wrap, p_wrap, ncol = 1, nrow = 3)
# P_wrap3
# 
# total_antho_combined <- ggarrange(p_key, p_ordi_totalantho, ncol = 1, nrow = 2, common.legend = TRUE, legend = "top")
# total_antho_combined
# 
# final_antho <- ggarrange(p_key, P_wrap3, ncol = 2, nrow = 1, align = "h", widths = c(2,1))
# final_antho

final_key_sp <- grid.arrange(legend_grob_test, p_key, ncol = 2, widths = c(1,5))
final_key_sp

ggsave(final_key_sp, file=filepath2, width = 40, height = 20, units = "cm", useDingbats=F)

# final_antho_con <- ggarrange(p_ordi_totalantho, anthos, ncol = 2, nrow = 1, align = "h", widths = c(2,1))
# final_antho_con

final_anti_sp_con <- grid.arrange(legend_grob_test, anti_combo, ncol = 2, widths = c(1,1))
final_anti_sp_con

ggsave(final_anti_sp_con, file=filepath3, width = 40, height = 20, units = "cm", useDingbats=F)


stop("Not an error - code done, anthocyanin plots generated")



###-----------End of code---------------------
