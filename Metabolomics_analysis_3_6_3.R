#Code written by M. R. Aggerbeck. 
#mrag@envs.au.dk
#Creative commons, 2023. 

###----------------PACKAGES--------------

library("tidyverse"); packageVersion("tidyverse")
library("ggpubr"); packageVersion("ggpubr")
library("data.table"); packageVersion("data.table")
library("igraph"); packageVersion("igraph")
library("plotrix"); packageVersion("plotrix")
library("VennDiagram"); packageVersion("VennDiagram")
library("phyloseq"); packageVersion("phyloseq")
library("devtools"); packageVersion("devtools")
library("biomformat"); packageVersion("biomformat")
library("vegan"); packageVersion("vegan")
library("DESeq2"); packageVersion("DESeq2")
library("mvabund"); packageVersion("mvabund")
library("metacoder"); packageVersion("metacoder")
library("taxa"); packageVersion("taxa")
library("adespatial"); packageVersion("adespatial")
library("microbiome"); packageVersion("microbiome")
library("pairwiseAdonis"); packageVersion("pairwiseAdonis")
library("microbiomeSeq"); packageVersion("microbiomeSeq")

#Save packages for easy .txt overview/snippet for paper m&m section
packages <- c(
  "tidyverse", "ggpubr", "data.table", "igraph", "plotrix",
  "VennDiagram", "phyloseq", "devtools", "biomformat", "vegan",
  "DESeq2", "mvabund", "metacoder", "taxa", "adespatial",
  "microbiome", "pairwiseAdonis", "microbiomeSeq"
)

packages_used <- data.frame(
  Package = character(),
  Version = character(),
  stringsAsFactors = FALSE
)

for (pkg in packages) {
  library(pkg, character.only = TRUE)
  version <- packageVersion(pkg)
  packages_used <- rbind(packages_used, data.frame(Package = pkg, Version = version))
}

packages_used
write.table(packages_used, "packages_used.txt", sep = "\t", quote = FALSE, row.names = FALSE)

###-----------Functions----------------

#----------ANOVA + pvalues------------------------------
plot_anova_diversity_pval <- function(physeq, method, grouping_column, bonf=TRUE, pValueCutoff=0.05)
{

  
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  abund_table <- otu_table(physeq)
  meta_table <- sample_data(physeq)
  
  #get diversity measure using selected methods
  div.df <- alpha_div(physeq,method)
  
  #=add grouping information to alpha diversity measures
  df<-data.frame(div.df,(meta_table[,grouping_column])[as.character(div.df$sample),])
  
  #perform anova of diversity measure between groups
  anova_res <- perform_anova(df,meta_table,grouping_column,pValueCutoff)
  df_pw <- anova_res$df_pw #get pairwise p-values
  
  #Apply Bonferroni correction
  if(bonf){  
    
    numberofsites <- length(unique(c(as.vector(levels(df_pw$from)),as.vector(levels(df_pw$to)))))[1]
    bonf.cor <- as.numeric(as.matrix(df_pw$p))* numberofsites
    temp <- as.factor(bonf.cor)
    df_pw$p <- temp
    }

  

  
  #Draw the boxplots
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
  p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  p<-p+theme_bw()
  p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Alpha Diversity Measure")+xlab("")
  p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("")
  

    
  #This loop will generate the lines and signficances
  if(!is.null(df_pw)){ #this only happens when we have significant pairwise anova results
    for(i in 1:dim(df_pw)[1]){
      p<-p+geom_path(inherit.aes=F,aes(x,y),
                     data = data.frame(
                       x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),
                             which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), 
                       y = c(as.numeric(as.character(df_pw[i,"y"])),
                             as.numeric(as.character(df_pw[i,"y"]))), 
                       measure=c(as.character(df_pw[i,"measure"]),
                                 as.character(df_pw[i,"measure"]))), 
                     color="black",lineend = "butt",
                     arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),
                     data=data.frame(
                       x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+
                            which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,
                       y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),
                       label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),
                                              breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),
                                              label=c("***", "**", "*", "")))))
    }
  }
  newlist <- list(p, df_pw)
  return(newlist)
}




###---------START OF CODE---------------

#Set working directory
setwd("O:/Tech_ENVS-EMBI-Afdelingsdrev/Marie/PhD/DATA/RESISTOME_QIIME/Analysis/Metabolomics")
#Set user directory
uzdir <- "O:/Tech_ENVS-EMBI-Afdelingsdrev/Marie/PhD/DATA/RESISTOME_QIIME/Analysis/Metabolomics"

# ###------------Find duplicates in preliminary dataset-------------------
# 
# 
# dup <- duplicated(NTS_prelim[,1])
# NTS <- NTS_prelim[,1]
# 
# is_dup <- NTS[duplicated(NTS)]
# 
# is_dup
# dim(is_dup)
# 
# duplicates <- unique(is_dup)                     
# write.csv2(duplicates, "duplicates.csv")
# 
# 
# #/
# 


###----------Import cleaned, annotated dataset from .csv---------------------
#Clean headers etc in txt editor/excel - remove regexes, special characters etc. 

# NTS <- read.csv2(file = "NTS_cleaned.csv")
RAW <- read.csv2(file = "Dataset_R_classifications_preliminary.csv", header = TRUE, sep=",")

# n<-dim(RAW)[1]
# RAW<-RAW[1:(n-10),]
dim(RAW)

head(RAW)
class(rownames(RAW))
colnames(RAW)[3]
colnames(RAW)

# RAW_sign

#Rownames will be unique identifiers for each compound - must be unique value! Use compound ID? 

rowname <- as.data.frame(RAW[1])
rowname <- as.character(as.matrix(RAW[1]))
rowname[is.na(rowname)] <- "unidentified"

head(rowname)
rowname <- paste("Compound", rowname, sep = "_")
rowname <- make.unique(rowname)

rawmetadata <- RAW

# for (i in seq_along(RAW[3, ])) {
#   if (RAW[3, i] == "" | is.na(RAW[3, i])) {
#     RAW[4, i] <- "Unknown Structure"
#   }
# }

colnames(RAW)
rownames(RAW) = rowname
rownames(RAW)


#Make compound table (OTU table equivalent - area/counts) 
#Find rows with areas. 
COMP_group <- as.data.frame(RAW[143:145])
COMP_group <- apply(COMP_group, 2, as.numeric)
head(COMP_group)
COMP_group[is.na(COMP_group)] = 0
colnames(COMP_group)
colnames(COMP_group) = c("Group_Area_A", "Group_Area_D", "Group_Area_S")

rownames(COMP_group) = rowname

head(COMP_group)
dim(COMP_group)

COMP <- as.data.frame(RAW[46:142])
rowname <- rownames(COMP)
COMP <- apply(COMP, 2, as.numeric)
head(COMP)
COMP[is.na(COMP)] = 0
colnames(COMP)
rownames(COMP) = rowname
rownames(COMP)
#Remove excess A_T2 samples
COMP <- COMP[,-c(7:12)]
head(COMP)

#Create preliminary metadat file from colname info
temp <- colnames(COMP) 
class(temp)
temp <- gsub("\\.", "å", temp)
temp <- gsub("åå", "å", temp)
temp <- data.frame(do.call(rbind, strsplit(temp, "å")))
temp

temp_meta <- as.character(temp$X2)
temp_meta <- gsub("A_", "Sample_A_", temp_meta)
temp_meta <- gsub("D_", "Sample_D_", temp_meta)
temp_meta <- gsub("S_", "Sample_S_", temp_meta)
temp_meta <- gsub("PB_Sample_", "PB_", temp_meta)
temp_meta <- gsub("QC", "QC_", temp_meta)
# temp_meta <- gsub("PB_", "PBå", temp_meta)
# temp_meta <- gsub("_T", "øT", temp_meta)
temp_meta <- strsplit(temp_meta, "_")
# temp_meta <- data.frame(do.call(rbind, strsplit(temp_meta, "_")))
temp_meta
max_elements <- max(sapply(temp_meta, length))
max_elements
split_matrix <- matrix("", nrow = length(temp_meta), ncol = max_elements)
for (i in 1:length(temp_meta)) {
  split_matrix[i, 1:length(temp_meta[[i]])] <- temp_meta[[i]]
}
result_df <- as.data.frame(split_matrix, stringsAsFactors = FALSE)
result_df
result_df <- cbind(temp$X2, result_df)
colnames(result_df) <- cbind("Sample_ID","Sample_Type", "Site", "Timepoint", "Bio_replicate", "Tech_replicate")
rownames(result_df) <- temp$X2
result_df 


colnames(COMP) = temp$X2
rownames(COMP) = rowname

head(COMP)
dim(COMP)

#Make drug class table (Tax table equivalent - "taxonomy" of drug classes - Therapeutics;Illegal drugs;Metabolite etc)
#Remember to pick the proper columns
CLA <- RAW[3:6]
rownames(CLA) = rowname
head(CLA)
dim(CLA)
colnames(CLA)
#Create rownames for metadata table

temp <- as.data.frame(colnames(COMP))
temp
#  
# #Make changes to sample IDs if too long/weird/etc. 
# temp$`colnames(COMP)` <- gsub('Area..','', temp$`colnames(COMP)`)
# temp$`colnames(COMP)` <- gsub('MRA.','', temp$`colnames(COMP)`)
# temp$`colnames(COMP)` <- gsub('WWTP.','', temp$`colnames(COMP)`)
# temp$`colnames(COMP)` <- gsub('raw..','', temp$`colnames(COMP)`, fixed=TRUE)
# temp #just double checking everything is in order
# temp$`colnames(COMP)` <- gsub('.PB','_PB', temp$`colnames(COMP)`, fixed=TRUE)
# temp$`colnames(COMP)` <- gsub('.QC','_QC', temp$`colnames(COMP)`, fixed=TRUE)
# temp #just double checking again
# temp$`colnames(COMP)` <- gsub('\\..*','', temp$`colnames(COMP)`)
# temp #triple check? 
# 
# #Swap names
# rownames(temp) = temp$`colnames(COMP)`
# colnames(COMP) =rownames(temp)
# # rownames(temp) = colnames(COMP)

#Make small metadata table for initial data exploration.

# site_abbr = c("A","D","S")
# site = c("Avedøre","Damhusåen","Skævinge")
# type = c("Sample","Sample","Sample","Procedural_Blank","Sample",
#          "Sample","Sample","Procedural_Blank","QC","QC","QC", "Sample","Sample","Sample","Procedural_Blank")
# PE = c(350000,350000,350000,350000,400000,400000,400000,400000,NA,NA,NA,12000,12000,12000,12000)
# temp2 <- cbind(temp,site, site_abbr)
# colnames(temp2) = c("Sample_ID","Sample_Site","Site") #,"Sample_Type","PE")
# temp2

result_df
#Add more columns if needed.
result_df$Site_Timepoint <- paste(result_df$Site, result_df$Timepoint, sep = "_")

metadata <- sample_data(result_df)
rownames(metadata) = metadata$Sample_ID
head(metadata)

COMP_phy = otu_table(COMP, taxa_are_rows = TRUE)
CLA_phy = tax_table(as.matrix(CLA))

any(is.na(COMP))
any(is.na(CLA))
any(is.na(metadata))


physeq = phyloseq(COMP_phy, CLA_phy)
physeq
physeq <- merge_phyloseq(physeq, metadata)
physeq

taxa_names(physeq) = rowname
head(taxa_names(physeq))

sample_names(physeq)

#For grouped areas: 
COMP_gr_phy = otu_table(COMP_group, taxa_are_rows = TRUE)

physeq_gr = phyloseq(COMP_gr_phy, CLA_phy)
physeq_gr
physeq_gr <- merge_phyloseq(physeq_gr, metadata)
physeq_gr

taxa_names(physeq_gr) = rowname
head(taxa_names(physeq_gr))

sample_names(physeq_gr)
head(tax_table(physeq_gr))
head(otu_table(physeq_gr))
head(metadata)

# #Code for removal of sample, in case of contamination
# physeq = subset_samples(physeq, Sample_ID != "All_samples" & Sample_ID != "All_samples.1")

#backup object in case of accidental erasure/change of primary physeq object
backup <- physeq
backup_gr <- physeq_gr
#Use this line to restore physeq - remember to run everything below this line again, to make sure nothing is messed up.
# physeq <- backup

# #Use if necessary:
# #Clean up character vectors - remove underscores, add hyphens, slashes etc. 
# 
# tax_table(physeq)[,1] <- sub("Textile/Chemicals/Auxiliary Dyes", 
# "Textile Chemicals/Auxiliary Dyes", tax_table(physeq)[,1])
# unique(tax_table(physeq)[,1])
# tax_table(physeq)[,2] <- sub("Fluorescent dye, possibly contaminant? ", "Fluorescent dye", tax_table(physeq)[,2])
# unique(tax_table(physeq)[,2])
# tax_table(physeq)[,2]

# physeq_noPB = subset_samples(physeq, Sample_ID != "Area_A_PB" & Sample_ID != "Area_D_PB" & Sample_ID != "Area_S_PB" )
# physeq_noPB
# 
# physeq_noQC = subset_samples(physeq_noPB, Sample_ID != "Area_Pool_QC_1" & 
#                                Sample_ID != "Area_Pool_QC_2" & Sample_ID != "Area_Pool_QC_3" )
# physeq_noQC
# 
# drugs = subset_taxa(physeq_noQC, Class=="Therapeutics/Drugs")
# drugs

No_PBQC <- subset_samples(physeq, Sample_Type == "Sample")
No_PBQC

No_PBQC_gr <- subset_samples(physeq_gr, Sample_Type == "Sample")
No_PBQC_gr

No_PB <- subset_samples(physeq, Sample_Type != "PB")
No_PB

unique(tax_table(physeq)[,2])

#Subset without unknowns
known <- subset_taxa(No_PBQC, !Class=="Unknown Structure")
# known <- subset_taxa(known, !Class=="Unknown Natural Product")
known
unique(tax_table(known)[,2])

IDonly <- subset_taxa(known, !Class=="Unknown Natural Product")
IDonly
unique(tax_table(IDonly)[,2])

# 
# #no blanks
# known_nopb <- subset_taxa(physeq_noPB, !Class=="Unknown structure")
# known_nopb
# unique(tax_table(known_nopb)[,2])
# 
# #samples only
# known_so <- subset_taxa(physeq_noQC, !Class=="Unknown structure")
# known_so
# unique(tax_table(known_so)[,2])
# # 
# #subset without unknown metabolites
# knownmet <- subset_taxa(known, !Class=="Unidentified Natural Product")
# knownmet
# knownmet_nopb <- subset_taxa(known_nopb, !Class=="Unidentified Natural Product")
# knownmet_nopb
# knownmet_so <- subset_taxa(known_so, !Class=="Unidentified Natural Product")
# knownmet_so
# 
# #test
# test_sunscr = subset_taxa(knownmet_so, Class=="Personal Care Products")
# tax_table(test_sunscr)

# #Write table with corrected names for substituting into the paper table
# 
# write.table(as.data.frame(tax_table(physeq)), file='tax_table.tsv', quote=FALSE, sep='\t')

###-------log-transform data

#physeq_noQC_log10 <- microbiome::transform(physeq_noQC, "log10p")

###-------Relative abundance

# physeq  = transform_sample_counts(physeq, function(x) x / sum(x) )

###--------Treatment Plant Relative Abundance (TPRA) ---------------

TPRA <- known_so

# test = prune_taxa(names(sort(taxa_sums(known_so), TRUE))[1:10], known_so) 
# test
# sample_data(test)
# otu_table(test)
# 
# test2 <- otu_table(test)
# otu_table(test)
# test2
# 
# for(n in 1:nsamples(test)){
#   otu_table(test)[,n] <- otu_table(test)[,n]/sample_data(test)$PE [n]
# }

for(n in 1:nsamples(TPRA)){
  otu_table(TPRA)[,n] <- otu_table(TPRA)[,n]/sample_data(TPRA)$PE [n]
}

###--------Merge samples into single group (viz only)------

temp = prune_taxa(taxa_sums(No_PBQC) > 0, No_PBQC)
temp2 <- as.data.frame(unique(sample_data(temp)[,7]))
groups_phy <- as.vector(as.character(temp2$Site_Timepoint))
sample_data(temp)$group_phy <- get_variable(temp, "Site_Timepoint") %in% groups_phy

merged_phy = merge_samples(temp, "Site_Timepoint")

merged_phy
sample_names(temp)
sample_names(merged_phy)

sample_data(merged_phy)
sample_data(merged_phy)$Site_Timepoint = sample_names(merged_phy)
merged_phy
physeq

tmp <- sample_data(merged_phy)

class(tmp)

tmp <- as.data.frame(tmp)
tmp
tmp2 <- tmp$Site_Timepoint
tmp2 <- str_split(tmp2, "_", simplify = TRUE); tmp2
tmp$Site <- tmp2[,1]
tmp$Timepoint <- tmp2[,2]
tmp

write.table(tmp, file='merged_metadata.tsv', quote=FALSE, sep='\t')

sample_data(merged_phy) <- tmp

known_merged <- subset_taxa(merged_phy, !Class=="Unknown Structure")
known_merged <- subset_taxa(known_merged, !Class=="Unknown")
known_merged

# ###------------Core compounds----------------
# NB! Does not work with this dataset - maybe set all below a certain threshold to 0? 

# 
# core.taxa.standard <- core_members(known_so, detection = 5/100, prevalence = 90/100)
# core.taxa.standard
# class(core.taxa.standard)
# write.table(core.taxa.standard, file='core_compounds.tsv', quote=FALSE, sep='\t')
# 
# Core_subset <- subset(otu_table(known_so), rownames(otu_table(known_so)) %in% core.taxa.standard)
# Core_subset
# new_physeq <- merge_phyloseq(my_subset, tax_table(physeq), sample_data(physeq), ...)

###--------tests-------------

# #----Richness - maybe?----
# 
# temp <- physeq_noQC
# head(otu_table(temp))
# otu_table(temp) = round(otu_table(temp))
# 
# plot_richness(temp, x = 'Site', measures=c("Observed", "Shannon", "Simpson"))+ 
#   theme_bw()+
#   theme(legend.title = element_blank())+
#   theme(legend.position = "none")+
#   theme(strip.background = element_rect(fill="white" ))
# 
# p_an <-plot_anova_diversity_pval(physeq_noQC, method = c("richness","shannon", "simpson"),
#                                  grouping_column ="Site", bonf=TRUE, pValueCutoff=0.05)
# 
# print(p_an)
# p <- p_an[[1]]
# p
# ggsave(p, file="RAW_noQC_richness.pdf",
#        width = 35, height = 20, units = "cm")

#----Barplot----

# Colour palette: 

set_class <- c("#8de985", "#412290","#019d4c","#9a0065","#01b6a8",
               "#ffacfd","#ae7400","#0195fb","#9a5184")
set_class10 <- c("#8de985", "#412290","#019d4c","#9a0065","#01b6a8",
               "#ffacfd","#ae7400","#0195fb","#9a5184", "black")

set_class11 <-  c("#8de985", "#412290","#019d4c","#9a0065","#01b6a8",
                  "#ffacfd","#ae7400","#0195fb","#9a5184", "white", "black")

# temp <- tax_glom(pseq.rel, taxrank = "Class")


p <- plot_bar(physeq_T6, fill="Class", title="Compounds found at all sites (except unknown structures)")+
  theme_bw() +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")+
  scale_fill_manual(values=set_class11)+
  scale_color_manual(values=set_class11)+
  labs(y= "Abundance", x = "Samples")
print(p)
# p2 <- p + facet_wrap(~Site+Timepoint, shrink = TRUE, drop = TRUE, 3)
p2 <- p + facet_grid(~Site+Timepoint, shrink = TRUE, drop=TRUE, scale="free")
p2
ggsave(p, file="RAW_noQC_barplot_Class.pdf", 
       width = 21, height = 14, units = "cm")

pseq.rel <- microbiome::transform(physeq_T6, "compositional")

p <- plot_bar(pseq.rel, fill="Class", title="Compounds found in samples, except unknown structures")+
  theme_bw() +
  labs(y= "Relative Abundance", x = "Samples")+
  scale_fill_manual(values=set_class11)+
  scale_color_manual(values=set_class11)+
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p)
ggsave(p, file="RAW_barplot_T6_Class_RelAb.pdf", 
       width = 21, height = 14, units = "cm")
p2 <- p + facet_grid(Site ~ Timepoint, shrink = TRUE, drop=TRUE, scales="free")+
  theme(axis.text.x = element_blank())
p2
ggsave(p2, file="RAW_barplot_T6_Class_RelAb_facet.pdf", 
       width = 21, height = 14, units = "cm")


#---ordination----

#Choosing bray-curtis, as dataset contains null values

# pseq.rel <- microbiome::transform(physeq_noPB, "compositional")

temp = prune_taxa(taxa_sums(No_PBQC) > 0, No_PBQC)

phy.ord <- ordinate(physeq_T6, "PCoA", "bray", pvalue.cutoff = 0.05)

p1 = plot_ordination(physeq_T6, phy.ord, type="sites", color="Site", shape="Timepoint", 
                     title="PCoA, Bray-Curtis Dissimilarity \n\n Timepoints")+
  geom_point(size=5)+
  scale_shape_manual(values=c(17, 15, 18, 19, 10, 5))+
  theme_bw()+
  # theme(legend.title = element_blank())+
  theme(#legend.position="bottom",
    strip.background = element_rect(fill="white" ))
print(p1)
p1_wrap <- p1 + facet_wrap(~Timepoint)
p1_wrap

ggsave(p1, file="RAW_Ordination_Samples_T6.pdf", 
       width = 30, height = 15, units = "cm", useDingbats=FALSE)


p2 = plot_ordination(physeq_T6, phy.ord, type="species", color="Class", 
                     title="\n\nB, Individual compounds")+
  geom_point(size=2, shape=16)+
  theme_bw()+
  scale_fill_manual(values=set_class10)+
  scale_color_manual(values=set_class10)+
  # theme(legend.title = element_blank())+
  facet_wrap(~Class, 3)+
  theme(#legend.position="bottom",
    strip.background = element_rect(fill="white" ))
print(p2)

# p1 <- p1 + scale_color_discrete(labels = c("Compound", "A", "D", "QC", "S")) +
#   scale_size_discrete(labels = c("Sample", "Compound")) +
#   scale_shape_discrete(labels = c("Compound", "QC", "Sample"))
# 
# print(p1)

# # Merging legends
# legend_1 <- get_legend(p1)
# legend_2 <- get_legend(p2)
# 
# legends <- ggarrange(legend_1, legend_2, nrow=2)
# 
# # Combining plots
# rm_legend <- function(p){p + theme(legend.position = "none")}
# plots <- ggarrange(rm_legend(p1), rm_legend(p2), nrow = 1, align = "v")
# 
# # plots + merged legends
# p3 <- ggarrange(plots, legends, widths = c(0.75, 0.25))

p3 <- ggarrange(p1,p2,#legend_1, legend_2,
               widths = 1, heights = 1.1,
               # labels = "PCoA, Bray-Curtis Dissimilarity",
               #label.x = 0.8, label.y = .90, hjust = 0, vjust = 0,
               ncol = 2, nrow = 1
               # common.legend = TRUE, legend = "right"
               )
p3

ggsave(p3, file="RAW_Ordination_SamplesvsQC_comp.pdf", 
       width = 30, height = 15, units = "cm", useDingbats=FALSE)

phy.ord <- ordinate(known_nopb, "PCoA", "bray", pvalue.cutoff = 0.05)
p1 = plot_ordination(known_nopb, phy.ord, type="Biplot", color="Sample_type", shape="Sample_Type", 
                     title="PCoA, bray-curtis dissimilarity, -unknowns")+
  geom_point(size=2)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))
print(p1)

ggsave(p1, file="RAW_Ordination_SamplesvsQC_col.pdf", 
       width = 25, height = 20, units = "cm", useDingbats=FALSE)


p2 <- p1 + facet_wrap(~Class, 3)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))

# p2 + geom_polygon(aes(fill=Sample_Type)) + geom_point(size=5)

print(p2)

ggsave(p2, file="RAW_Ordination_SamplesvsQC_col_wrapped.pdf", 
       width = 25, height = 20, units = "cm", useDingbats=FALSE)


phy.ord <- ordinate(knownmet_so, "PCoA", "bray", pvalue.cutoff = 0.05)
p1 = plot_ordination(knownmet_so, phy.ord, type="Split", color="Class", shape="Site", 
                     title="PCoA, bray-curtis dissimilarity")+
  geom_point(size=2)+
  scale_shape_manual(values=c(16,15,17,8))+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))
print(p1)

ggsave(p1, file="RAW_Ordination_PCoA_bray_Class.pdf", 
       width = 35, height = 20, units = "cm", useDingbats=FALSE)


p2 <- p1 + facet_wrap(~Class, 3)+
  theme_bw()+
  theme(legend.title = element_blank())+
  theme(
    strip.background = element_rect(fill="white" ))

# p2 + geom_polygon(aes(fill=Sample_Type)) + geom_point(size=5)

print(p2)

ggsave(p2, file="RAW_Ordination_PCoA_bray_Class_fw.pdf", 
       width = 35, height = 20, units = "cm", useDingbats=FALSE)


#---Heatmap

tmp <- known

tmp <- prune_taxa(names(sort(taxa_sums(tmp),TRUE)[1:50]), tmp)
p <- plot_heatmap(tmp, method="NMDS", distance="bray", sample.label="Site", taxa.label = "Name", 
                  sample.order = "Sample_Type", title="Heatmap, NMDS + Bray-Curtis") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank()) +
  labs(fill = "Area")
print(p)
ggsave(p, file="RAW_noQC_Heatmap.pdf", width = 35, height = 20, units = "cm")

tmp <- prune_taxa(names(sort(taxa_sums(tmp),TRUE)[1:50]), tmp)
p <- plot_heatmap(tmp, method="NMDS", distance="bray", sample.label="Site", 
                  sample.order = "Sample_Type", title="Heatmap, NMDS + Bray-Curtis, 50 most abundant compounds") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank()) +
  labs(fill = "Area")
print(p)
ggsave(p, file="RAW_noQC_Heatmap_50mostab.pdf", width = 35, height = 20, units = "cm")


###--------PERMANOVA + post-hoc pairwise---------------------------

# Insert object
pseq <- physeq_S
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- microbiome::abundances(pseq)
meta <- microbiome::meta(pseq)


permanova_S <- vegan::adonis(t(otu)~ Timepoint,
                           data = meta, permutations=999, method = "bray")


# permanova$aov.tab

post_hoc_permanova_S <- pairwise.adonis(t(otu), meta$Timepoint, sim.function = "vegdist",
                                        sim.method = "bray", p.adjust.m = "fdr", reduce = NULL,
                                        perm = 999)
#Anova Results
#Site A
permanova_A$aov.tab
#Site D
permanova_D$aov.tab
#Site S
permanova_S$aov.tab


# post_hoc_permanova
#Site A
post_hoc_permanova_A
#Site D
post_hoc_permanova_D
#Site S
post_hoc_permanova_S


###--------DESeq2------------

#Be aware that Deseq2 has a min and a max value. min=0, max= 2.9e9, 
#Which is occasionally exceeded when using chemical data (area). 

#DESeq also gets angry if you log transform. Scale down all counts by a factor 10 or 100, to bypass this error.

head(otu_table(physeq_noQC))
test_phy <- (otu_table(physeq_noQC))/100
head(test_phy)

#DESeq2

temp <- physeq_noQC
otu_table(temp) = test_phy
head(otu_table(temp))

diagdds = phyloseq_to_deseq2(temp, ~Site)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
print(diagdds)

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(temp)[rownames(sigtab), ], "matrix"))
head(sigtab)

# theme_set(theme_bw())

# Class order
x = tapply(sigtab$log2FoldChange, sigtab$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab$Class = factor(as.character(sigtab$Class), levels=names(x))
# Subclass order
x = tapply(sigtab$log2FoldChange, sigtab$Subclass_I, function(x) max(x))
x = sort(x, TRUE)
sigtab$Subclass_I = factor(as.character(sigtab$Subclass_I), levels=names(x))
ggplot(sigtab, aes(x=Subclass_I, y=log2FoldChange, color=Class)) + geom_point(size=2) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# OTU <- otu_table(physeq_noQC, taxa_are_rows = TRUE)
# OTU <- t(OTU)
# #Merge new phyloseq object
# TAX = as.matrix(tax_table(physeq_noQC))
# met <- as.data.frame(sample_data(physeq_noQC))
# pseq1 <- phyloseq(OTU, TAX)
# physeq2<-merge_phyloseq(pseq1, met)
# rm(OTU, TAX, met, pseq1)
# 
# minTotRelAbun = 0.1
# x = taxa_sums(physeq2)
# keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
# phy3 = prune_taxa(keepTaxa, physeq2)
# phy3 = physeq2
# 
# deseq_sig <- differential_abundance(physeq_noQC, grouping_column = "Site",
#                                     pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = "RAW_site_DeSeq2")


###---------Heat tree-------------------


# temp <- drugs
temp <- physeq_T6
# temp <- physeq_noQC

tax_table(temp) <- tax_table(temp)[,2:4]
head(tax_table(temp))

# temptax <- as.data.frame(tax_table(temp))
# head(temptax)
# temptax$Class <- gsub("^$", "Unclassified_compound", temptax$Class)
# temptax$Subclass <- gsub("^$", NA, temptax$Subclass)
# temptax$SubclassII <- gsub("^$", NA, temptax$SubclassII)
# temptax = tax_table(as.matrix(temptax))
# tax_table(temp) <- temptax
# temp <- tax_glom(temp, taxrank = "Subclass")


#parse phyloseq object temp

obj_all <- parse_phyloseq(temp)

obj <- obj_all
tissuegroup <- obj$data$sample_data$Site

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")


heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "fruchterman-reingold",
                 initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'RAW_heattree_knownmet_T6.pdf')


# x = parse_tax_data(hmp_otus, class_cols = "lineage", class_sep = ";",
#                    class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
#                    class_regex = "^(.+)__(.+)$")
# 
# x
# # Get per-taxon counts
# x$data$tax_table <- calc_taxon_abund(x, data = "tax_data", cols = hmp_samples$sample_id)
# 
# # Calculate difference between groups
# x$data$diff_table <- calc_diff_abund_deseq2(x, data = "tax_table",
#                                             cols = hmp_samples$sample_id,
#                                             groups = hmp_samples$body_site)

# Remove taxa with only small differences 
per_taxon_fold_changes <- obs(obj, data = 'diff_table', value = 'log2_median_ratio')
per_taxon_fold_changes
per_taxon_max_change <- unlist(lapply(per_taxon_fold_changes, function(tax_changes) max(abs(tax_changes))))
per_taxon_max_change
x <- filter_taxa(obj, per_taxon_max_change > 3, supertaxa = TRUE, reassign_obs = c(diff_table = FALSE))
x

taxa_list_classification <- as.data.frame(x$data$tax_data)
taxa_list_classification

write.table(taxa_list_classification, file='taxa_list_classification.tsv', quote=FALSE, sep='\t')

heat_tree_matrix(x,
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
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'RAW_heattree_LfC1.pdf')


# #Test p-values
# # obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
# test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.1] <- 0 )
# obj$data$diff_table <- test
# 
# heat_tree_matrix(obj,
#                  data = "diff_table",
#                  node_size = n_obs,
#                  node_label = taxon_names,
#                  node_color = log2_median_ratio,
#                  node_color_range = diverging_palette(),
#                  node_color_trans = "linear",
#                  node_color_interval = c(-3, 3),
#                  edge_color_interval = c(-3, 3),
#                  layout = "davidson-harel",
#                  initial_layout = "reingold-tilford",
#                  node_size_axis_label = "Number of compounds",
#                  node_color_axis_label = "Log2 ratio median proportions",
#                  output_file = 'RAW_heattree_pval.pdf')


###---------Venn Diagram----------------------------

Venn <- read.csv2(file = "venn.csv") 
colnames(Venn)

overlap <- calculate.overlap(x = list('A' = set1, 'D' = set2, 'S' = set3))
area_overlap <- sapply(overlap, length)

## areas of each animal
area1 <- as.numeric(nrow(subset(Venn, A == 1)))
## [1] 135
area2 <- as.numeric(nrow(subset(Venn, D == 1)))
## [1] 99
area3 <- as.numeric(nrow(subset(Venn, S == 1)))
## [1] 92


## areas of 2-group overlap

## A & D
n12 <- as.numeric(nrow(subset(Venn, A == 1 & D == 1)))
## [1] 75

## A & S
n13 <- as.numeric(nrow(subset(Venn, A == 1 & S == 1)))
## [1] 52

## D & S
n23 <- as.numeric(nrow(subset(Venn, D == 1 & S == 1)))
## [1] 53

## 3 group overlap: A, D & S
n123 <- as.numeric(nrow(subset(Venn, A == 1 & D == 1 & S == 1)))
## [1] 47

draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category =
                   c("A", "D", "S"), lty = "blank", 
                 col = rep("black", 3), fill = c("darkblue", "darkorange", "turquoise"))

###------------Sites--------------------

physeq_A = subset_samples(No_PBQC, Site == "A")
physeq_A = prune_taxa(taxa_sums(physeq_A) > 0, physeq_A)
physeq_A

physeq_D = subset_samples(No_PBQC, Site == "D")
physeq_D = prune_taxa(taxa_sums(physeq_D) > 0, physeq_D)
physeq_D

physeq_S = subset_samples(No_PBQC, Site == "S")
physeq_S = prune_taxa(taxa_sums(physeq_S) > 0, physeq_S)
physeq_S


###------------Timepoints---------------

physeq_T6 = subset_samples(No_PBQC, Timepoint == "T6")
physeq_T6 = prune_taxa(taxa_sums(physeq_T6) > 0, physeq_T6)
physeq_T6

###------------Sign only-----------------

RAW_sign <- RAW[RAW$Tag_sign =="S",]
RAW_sign

dim(RAW)
dim(RAW_sign)
head(RAW_sign)

rownames_sign <- as.character(rownames(RAW_sign))

physeq_sign <- prune_taxa(rownames_sign, known)
physeq_sign

timepoints <- sample_data(physeq_sign)$Timepoint
timepoints





###Merge first? 

temp = prune_taxa(taxa_sums(physeq_sign) > 0, physeq_sign)
temp2 <- as.data.frame(unique(sample_data(temp)[,7]))
groups_phy <- as.vector(as.character(temp2$Site_Timepoint))
sample_data(temp)$group_phy <- get_variable(temp, "Site_Timepoint") %in% groups_phy
merged_phy_sign = merge_samples(temp, "Site_Timepoint"); merged_phy_sign
sample_names(temp);sample_names(merged_phy_sign)

sample_data(merged_phy_sign)
sample_data(merged_phy_sign)$Site_Timepoint = sample_names(merged_phy_sign)
merged_phy_sign; physeq_sign
tmp <- sample_data(merged_phy_sign)
tmp <- as.data.frame(tmp); tmp
tmp2 <- tmp$Site_Timepoint
tmp2 <- str_split(tmp2, "_", simplify = TRUE); tmp2
tmp$Site <- tmp2[,1]
tmp$Timepoint <- tmp2[,2]
tmp
write.table(tmp, file='merged_sign_metadata.tsv', quote=FALSE, sep='\t')
sample_data(merged_phy_sign) <- tmp

# known_merged <- subset_taxa(merged_phy_sign, !Class=="Unknown Structure")
# known_merged <- subset_taxa(known_merged, !Class=="Unknown")
# known_merged

#Barplot

ps <- plot_bar(merged_phy_sign, fill="Class", title="Significant compounds in samples (except unknown structures)")+
  theme_bw() +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")+
  scale_fill_manual(values=set_class)+
  scale_color_manual(values=set_class)+
  # scale_x_discrete(labels = NULL) +
  labs(y= "Abundance", x = "Timepoint")
print(ps)
p2 <- ps + facet_wrap(~Site, nrow = 3, scales = "free")
# p2 <- ps + facet_wrap(~Site, shrink = TRUE, drop = TRUE, scale="free", 3)
p2 <- ps + facet_grid(~Timepoint + Site, shrink = TRUE, drop=TRUE, scales="free", space = "free")
p2
ggsave(ps, file="RAW_sign_noQC_barplot_Class.pdf", 
       width = 21, height = 14, units = "cm")

pseq.rel <- microbiome::transform(merged_phy_sign, "compositional")

ps_rel <- plot_bar(pseq.rel, fill="Class", title="Significant compounds in samples (except unknown structures)")+
  theme_bw() +
  labs(y= "Relative Abundance", x = "Samples")+
  scale_fill_manual(values=set_class)+
  scale_color_manual(values=set_class)+
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
print(ps_rel)
ps_rel <- ps_rel + facet_wrap(~Site, ncol = 3, scales = "free")
ggsave(ps_rel, file="RAW_sign_barplot_Class_RelAb.pdf", 
       width = 21, height = 14, units = "cm")


sample_data(physeq_sign)

###------Antibiotics-----

ps_anti <- subset_taxa(temp, Class == "Therapeutics/Drugs")
ps_anti

tax_table(ps_anti)
otu_table(ps_anti)


# pseq.rel <- microbiome::transform(ps_anti, "compositional")

ps_rel <- plot_bar(ps_anti, fill="Subclass", title="Antibiotic compounds")+
  theme_bw() +
  labs(y= "Relative Abundance", x = "Samples")+
  # scale_fill_manual(values=set_class)+
  # scale_color_manual(values=set_class)+
  geom_bar(aes(color=Subclass, fill=Subclass), stat="identity", position="stack")
print(ps_rel)
ps_rel + facet_wrap(~Timepoint, ncol = 3, scales = "free_x")
ggsave(ps_rel, file="RAW_sign_barplot_Class_RelAb.pdf", 
       width = 21, height = 14, units = "cm")


###--------CECs-------

# Subsetting for multiple categories using %in%

#Need raw values for ordisurf, not relative abundance
# temp <- microbiome::transform(known, "compositional")

ps_CECs <- subset_taxa(known, Class %in% c("Therapeutics/Drugs", "Industrial Chemicals", "Personal Care Products", "Pesticides"))
ps_CECs

tax_table(ps_CECs)[tax_table(ps_CECs) == "UNII:K0N34Z34O8"] <- "PPG-3 Butyl Ether"
tax_table(ps_CECs)[5,3] = "Polymer"; tax_table(ps_CECs)[5,4] = "PPG"

tax_table(ps_CECs)
otu_table(ps_CECs)
sample_data(ps_CECs)

CEC_names <- as.data.frame(tax_table(ps_CECs)[,1])
CEC_names

CEC_names$Name <- (str_replace_all(CEC_names$Name, "[^[:alnum:]]", "_"))
CEC_names$Name <- (str_replace_all(CEC_names$Name, "_+", "_"))

rownames(tax_table(ps_CECs))

taxa_names(ps_CECs) <- make.unique(CEC_names$Name)

plot_CECs <- plot_bar(ps_CECs, fill="Class", title="Compounds of Emerging Concern")+
  theme_bw() +
  labs(y= "Relative Abundance", x = "Samples")+
  scale_fill_manual(values=set_class)+
  scale_color_manual(values=set_class)+
  geom_bar(aes(color=Class, fill= Class), stat="identity", position="stack")
print(plot_CECs)
plot_CECs + facet_wrap(~Timepoint, ncol = 3, scales = "free_x")
ggsave(plot_CECs, file="CECs_barplot_Class_RelAb.pdf", 
       width = 21, height = 14, units = "cm")

###-------Test network---------------

set.seed(1)

ig <- make_network(ps_CECs, type="samples", max.dist=0.4)
ig_plot <- plot_network(ig, physeq=ps_CECs, type="samples", 
                        color="Site", shape="Timepoint", point_size=4, alpha=1,
                        label="value", hjust = 1.35, 
                        line_weight=0.5, line_alpha=0.4,
                        layout.method=layout.fruchterman.reingold, title=NULL)
ig_plot

ig2 <- make_network(ps_CECs, type="taxa", max.dist=0.4)
set.seed(1)
ig_plot1 <- plot_network(ig2, physeq=ps_CECs, type="taxa", 
                         color="Class", point_size=4, alpha=1,
                         label="value", hjust = 1.35, 
                         line_weight=0.5, line_alpha=0.4,
                         layout.method=layout.fruchterman.reingold, title=NULL)
ig_plot1

#----- separate site networks------------------

physeqCECs_A = subset_samples(ps_CECs, Site == "A")
physeqCECs_A = prune_taxa(taxa_sums(physeqCECs_A) > 0, physeqCECs_A)
physeqCECs_A

physeqCECs_D = subset_samples(ps_CECs, Site == "D")
physeqCECs_D = prune_taxa(taxa_sums(physeqCECs_D) > 0, physeqCECs_D)
physeqCECs_D

physeqCECs_S = subset_samples(ps_CECs, Site == "S")
physeqCECs_S = prune_taxa(taxa_sums(physeqCECs_S) > 0, physeqCECs_S)
physeqCECs_S


ig3 <- make_network(physeqCECs_A, type="taxa", max.dist=0.4)
set.seed(711L)
ig3_plot <- plot_network(ig3, physeq=physeqCECs_A, type="taxa", 
                         color="Class", point_size=4, alpha=1,
                         label="value", hjust = 1.35, 
                         line_weight=0.5, line_alpha=0.4,
                         layout.method=layout.fruchterman.reingold, title=NULL)
ig3_plot

ig4 <- make_network(physeqCECs_D, type="taxa", max.dist=0.4)
set.seed(711L)
ig4_plot <- plot_network(ig4, physeq=physeqCECs_D, type="taxa", 
                         color="Class", point_size=4, alpha=1,
                         label="value", hjust = 1.35, 
                         line_weight=0.5, line_alpha=0.4,
                         layout.method=layout.fruchterman.reingold, title=NULL)
ig4_plot

ig5 <- make_network(physeqCECs_S, type="taxa", max.dist=0.4)
set.seed(711L)
ig5_plot <- plot_network(ig5, physeq=physeqCECs_S, type="taxa", 
                         color="Class", point_size=4, alpha=1,
                         label="value", hjust = 1.35, 
                         line_weight=0.5, line_alpha=0.4,
                         layout.method=layout.fruchterman.reingold, title=NULL)
ig5_plot

ps_CECs
tax_table(ps_CECs)
otu_table(ps_CECs)

ps_CECs_ht <- ps_CECs
tax_table(ps_CECs_ht)[,1] <- "Compound"
tax_table(ps_CECs_ht)
tax_table(ps_CECs_ht)[,4] <- gsub("0", "", tax_table(ps_CECs_ht)[,4])
tax_table(ps_CECs_ht)[,4]

ps_CECs_ht_cl <- tax_glom(ps_CECs_ht, taxrank = "Subclass", NArm = FALSE)
ps_CECs_ht_cl
tax_table(ps_CECs_ht_cl)[,3]
taxa_names(ps_CECs_ht_cl) <- make.unique(tax_table(ps_CECs_ht_cl)[,3])
taxa_names(ps_CECs_ht_cl)

ps_CECs_ht_c <- tax_glom(ps_CECs_ht, taxrank = "Class", NArm = FALSE)
ps_CECs_ht_c
tax_table(ps_CECs_ht_c)[,2]
taxa_names(ps_CECs_ht_c) <- make.unique(tax_table(ps_CECs_ht_c)[,2])
taxa_names(ps_CECs_ht_c)


otu_1 <- t(as.data.frame(otu_table(ps_CECs_ht))[, sample_data(ps_CECs_ht)$Timepoint == "T6"])
otu_2 <- t(as.data.frame(otu_table(ps_CECs_ht_cl))[, sample_data(ps_CECs_ht_cl)$Timepoint == "T6"])
otu_3 <- t(as.data.frame(otu_table(ps_CECs_ht_c))[, sample_data(ps_CECs_ht_c)$Timepoint == "T6"])

head(otu_1)
head(otu_2)
head(otu_3)
otu_metabolites <- cbind(otu_3, otu_2, otu_1)
dim(otu_metabolites)
write.table(otu_metabolites, file='otu_metabolites_t6.tsv', quote=FALSE, sep='\t')
save(otu_metabolites, file = "otu_metabolites_t6.RData")



###---------Heat tree-------------------


# temp <- drugs
temp <- ps_CECs_ht
# temp <- physeq_noQC

# tax_table(temp) <- tax_table(temp)[,2:4]
# head(tax_table(temp))

temptax <- as.data.frame(tax_table(temp))
head(temptax)
# temptax$Class <- gsub("^$", "Unclassified_compound", temptax$Class)
# temptax$Subclass <- gsub("^$", NA, temptax$Subclass)
temptax$SubclassII <- gsub("0", "", temptax$SubclassII)
temptax = tax_table(as.matrix(temptax))
tax_table(temp) <- temptax
# temp <- tax_glom(temp, taxrank = "Subclass")


#parse phyloseq object temp

obj_all <- parse_phyloseq(temp)

obj <- obj_all
tissuegroup <- obj$data$sample_data$Site

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", 
                                      cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")


heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "fruchterman-reingold",
                 initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'CECs_heattree.pdf')



# Remove taxa with only small differences 
per_taxon_fold_changes <- obs(obj, data = 'diff_table', value = 'log2_median_ratio')
per_taxon_fold_changes
per_taxon_max_change <- unlist(lapply(per_taxon_fold_changes, function(tax_changes) max(abs(tax_changes))))
per_taxon_max_change
x <- filter_taxa(obj, per_taxon_max_change > 3, supertaxa = TRUE, reassign_obs = c(diff_table = FALSE))
x

taxa_list_classification <- as.data.frame(x$data$tax_data)
taxa_list_classification

write.table(taxa_list_classification, file='taxa_list_classification.tsv', quote=FALSE, sep='\t')

heat_tree_matrix(x,
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
                 node_size_axis_label = "Number of compounds",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'RAW_heattree_LfC1.pdf')


