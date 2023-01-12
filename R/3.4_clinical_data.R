library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(ggpubr)
library(gridExtra)
library(metadeconfoundR)
library(reshape2)

# ------------- DATA PREPROCESSING -----------------

# Read in data
emma_rds_file <- "ps_emma.rds"
ps_emma <- readRDS(emma_rds_file)

emma_metafile_baseline <- "220214_EMMA Datenbank-BL.csv"
meta_emma_bl <- read.csv(emma_metafile_baseline, sep = "\t", skip = 2)[-1,]
emma_metafile_2M <- "220214_EMMA Datenbank-2M.csv"
meta_emma_2M <- read.csv(emma_metafile_2M, sep = "\t", skip = 2)[-1,]
meta_emma_2M$ID[meta_emma_2M$ID == "01.2.S.030"] <- "01.1.S.030"
emma_metafile_4M <- "220214_EMMA Datenbank-4M.csv"
meta_emma_4M <- read.csv(emma_metafile_4M, sep = "\t", skip = 2)[-1,]
meta_emma_4M$ID[meta_emma_4M$ID == "01.2.S.030"] <- "01.1.S.030"

# Make Variables binary
levels(sample_data(ps_emma)$Sectio) <- c("Vaginal","C-section")
sample_data(ps_emma)$X2M_Ausschlag <- factor(sample_data(ps_emma)$X2M_Ausschlag)
levels(sample_data(ps_emma)$X2M_Ausschlag)[1] <- 888
levels(sample_data(ps_emma)$X2M_Ausschlag)[3:8] <- 1
sample_data(ps_emma)$X2M_dry_skin <- factor(sample_data(ps_emma)$X2M_dry_skin)
levels(sample_data(ps_emma)$X2M_dry_skin)[2:5] <- 1
sample_data(ps_emma)$X2M_Dx_eczem <- factor(sample_data(ps_emma)$X2M_Dx_eczem)
levels(sample_data(ps_emma)$X2M_Dx_eczem)[2:3] <- 1
sample_data(ps_emma)$X4M_Fieber <- factor(sample_data(ps_emma)$X4M_Fieber)
levels(sample_data(ps_emma)$X4M_Fieber)[c(1,4)] <- 0
levels(sample_data(ps_emma)$X4M_Fieber)[2:3] <- 1
sample_data(ps_emma)$X4M_Ausschlag <- factor(sample_data(ps_emma)$X4M_Ausschlag)
levels(sample_data(ps_emma)$X4M_Ausschlag)[2:3] <- 1
sample_data(ps_emma)$X4M_dry_skin <- factor(sample_data(ps_emma)$X4M_dry_skin)
levels(sample_data(ps_emma)$X4M_dry_skin)[2:5] <- 1
sample_data(ps_emma)$X4M_Dx_eczem <- factor(sample_data(ps_emma)$X4M_Dx_eczem)
levels(sample_data(ps_emma)$X4M_Dx_eczem)[2:4] <- 1
sample_data(ps_emma)$Ern_M_vege <- factor(sample_data(ps_emma)$Ern_M)
levels(sample_data(ps_emma)$Ern_M_vege)[c(1,4,5)] <- 0
levels(sample_data(ps_emma)$Ern_M_vege)[2:3] <- 1
sample_data(ps_emma)$Kompl_Geburt <- factor(sample_data(ps_emma)$Kompl_Geburt)
levels(sample_data(ps_emma)$Kompl_Geburt)[2:3]<-1
sample_data(ps_emma)$BG_M_A <- sample_data(ps_emma)$BG_M %in% c("1","2")
sample_data(ps_emma)$BG_M_B <- sample_data(ps_emma)$BG_M %in% c("2","3")
sample_data(ps_emma)$Geschw <- factor(sample_data(ps_emma)$Geschw)
levels(sample_data(ps_emma)$Geschw)[2:3] <- 1
sample_data(ps_emma)$Asthma_M <- factor(sample_data(ps_emma)$Asthma_M)
levels(sample_data(ps_emma)$Asthma_M)[c(1,3)] <- 0
sample_data(ps_emma)$Bronch_M <- factor(sample_data(ps_emma)$Bronch_M)
levels(sample_data(ps_emma)$Bronch_M)[c(1,3)] <- 0
sample_data(ps_emma)$Hay_M <- factor(sample_data(ps_emma)$Hay_M)
levels(sample_data(ps_emma)$Hay_M)[c(1,3)] <- 0
sample_data(ps_emma)$ND_M <- factor(sample_data(ps_emma)$ND_M)
levels(sample_data(ps_emma)$ND_M)[c(1,3)] <- 0
sample_data(ps_emma)$Ekz_M <- factor(sample_data(ps_emma)$Ekz_M)
levels(sample_data(ps_emma)$Ekz_M)[c(1,3)] <- 0
sample_data(ps_emma)$ImD_M <- factor(sample_data(ps_emma)$ImD_M)
levels(sample_data(ps_emma)$ImD_M)[c(2,3)] <- 1
sample_data(ps_emma)$All_M <- factor(sample_data(ps_emma)$All_M)
levels(sample_data(ps_emma)$All_M)[c(1,3)] <- 0
sample_data(ps_emma)$X2M_Fieber <- factor(sample_data(ps_emma)$X2M_Fieber)
levels(sample_data(ps_emma)$X2M_Fieber)[c(1,3)] <- 0
sample_data(ps_emma)$X2M_VitD <- factor(sample_data(ps_emma)$X2M_VitD)
levels(sample_data(ps_emma)$X2M_VitD)[1] <- 888
levels(sample_data(ps_emma)$X2M_VitD)[3:5] <- 1
sample_data(ps_emma)$X2M_NT <- factor(sample_data(ps_emma)$X2M_NT)
levels(sample_data(ps_emma)$X2M_NT)[c(1,3)] <- 0
sample_data(ps_emma)$X2M_rhi_nocold <- factor(sample_data(ps_emma)$X2M_rhi_nocold)
levels(sample_data(ps_emma)$X2M_rhi_nocold)[c(1,3)] <- 0
sample_data(ps_emma)$X1Impf_6fach <- factor(sample_data(ps_emma)$X1Impf_6fach)
levels(sample_data(ps_emma)$X1Impf_6fach)[1] <- 888
levels(sample_data(ps_emma)$X1Impf_6fach)[3:4] <- 1
sample_data(ps_emma)$X4M_VitD <- factor(sample_data(ps_emma)$X4M_VitD)
levels(sample_data(ps_emma)$X4M_VitD)[1] <- 888
levels(sample_data(ps_emma)$X4M_VitD)[3:5] <- 1
sample_data(ps_emma)$X4M_NT <- factor(sample_data(ps_emma)$X4M_NT)
levels(sample_data(ps_emma)$X4M_NT)[c(1,3)] <- 0
sample_data(ps_emma)$X2Impf_6fach <- factor(sample_data(ps_emma)$X2Impf_6fach)
levels(sample_data(ps_emma)$X2Impf_6fach)[2:3] <- 1

# Sum up variables
mother_immune_disorders<- c("Asthma_M","Bronch_M","Hay_M","ND_M","Ekz_M","ImD_M","All_M")
sample_data(ps_emma)$Immune_Disorder_M <- apply(sample_data(ps_emma)[,mother_immune_disorders],1,max,na.rm=T)

respiratory <- c("Erk","Bron","Grip","Pne","Kru","Oti","Ang","Herp","Konj","Hu","Hu_nachts","LN","Wh","rhi_nocold","eyes_itchy")
sample_data(ps_emma)$X2M_respiratory <- apply(sample_data(ps_emma)[,paste("X2M_",respiratory,sep="")],1,max,na.rm=T)
sample_data(ps_emma)$X4M_respiratory <- apply(sample_data(ps_emma)[,paste("X4M_",respiratory,sep="")],1,max,na.rm=T)

gastrointestinal <- c("GE","diar","HWI")
sample_data(ps_emma)$X2M_gastrointestinal <- apply(sample_data(ps_emma)[,paste("X2M_",gastrointestinal,sep="")],1,max,na.rm=T)
sample_data(ps_emma)$X4M_gastrointestinal <- apply(sample_data(ps_emma)[,paste("X4M_",gastrointestinal,sep="")],1,max,na.rm=T)

skin <- c("soor","Ausschlag","dry_skin","Dx_eczem")
sample_data(ps_emma)$X2M_skin <- apply(sample_data(ps_emma)[,paste("X2M_",append(skin,"Akne"),sep="")],1,max,na.rm=T)
sample_data(ps_emma)$X4M_skin <- apply(sample_data(ps_emma)[,paste("X4M_",skin,sep="")],1,max,na.rm=T)

vacc_2M <- c("X1Impf_6fach","X1Impf_Pneu","X1Impf_Rota")
vacc_4M <- append(vacc_2M,c("X1Impf_MenB","X2Impf_6fach","X2Impf_Pneu","X2Impf_Rota"))
sample_data(ps_emma)$X2M_vacc <- apply(sample_data(ps_emma)[,vacc_2M],1,max,na.rm=T)
sample_data(ps_emma)$X4M_vacc <- apply(sample_data(ps_emma)[,vacc_4M],1,max,na.rm=T)


# "888" to NA
sample_data(ps_emma)[sample_data(ps_emma) == "888"] <- NA

# create phyloseq objects
ps_rel <- transform_sample_counts(ps_emma,function(x) x/sum(x))

ps_visit1 <- prune_samples(get_variable(ps_emma, "visit") == "Visit 1" , ps_emma)
ps_visit2 <- prune_samples(get_variable(ps_emma, "visit") == "Visit 2" , ps_emma)
ps_visit3 <- prune_samples(get_variable(ps_emma, "visit") == "Visit 3" , ps_emma)

ps_v1_rel <- transform_sample_counts(ps_visit1, function(x) x/sum(x))
ps_v2_rel <- transform_sample_counts(ps_visit2, function(x) x/sum(x))
ps_v3_rel <- transform_sample_counts(ps_visit3, function(x) x/sum(x))

# ---------------- PERMANOVA -----------------

metadata <- as(sample_data(ps_rel),"data.frame")
features <- as.data.frame(otu_table(ps_rel))

metadata_v1 <- as(sample_data(ps_v1_rel),"data.frame")
features_v1 <- as.data.frame(otu_table(ps_v1_rel))

metadata_v2 <- as(sample_data(ps_v2_rel),"data.frame")
features_v2 <- as.data.frame(otu_table(ps_v2_rel))

metadata_v3 <- as(sample_data(ps_v3_rel),"data.frame")
features_v3 <- as.data.frame(otu_table(ps_v3_rel))

# ID "01.1.V.068" --> no metadata - only present in visit 1: 16011V068c001
metadata <- metadata[!(rownames(metadata)=="16011V068c001"),]
features <- features[!(rownames(features)=="16011V068c001"),]
metadata_v1 <- metadata_v1[!(rownames(metadata_v1)=="16011V068c001"),]
features_v1 <- features_v1[!(rownames(features_v1)=="16011V068c001"),]
# ID "01.2.S.089" --> no metadata for tot4M
metadata_v3 <- metadata_v3[!(rownames(metadata_v3)=="18012S089c090"),]
features_v3 <- features_v3[!(rownames(features_v3)=="18012S089c090"),]

set.seed(100)
perm_all <- adonis2(features ~ run+patient_id+visit+total_reads,data=metadata, permutations = 999, method = "bray", by = "margin")
write.table(perm_all,"3.4_PERMANOVA_all.csv",sep=",",col.names = NA)

set.seed(100)
perm_v1 <- adonis2(features_v1 ~ Sectio+sex+Ern_M_vege+ATB_SS+Bonding_KRS+MB_Kontakt,data=metadata_v1, permutations = 999, method = "bray", by = "margin")
write.table(perm_v1,"3.4_PERMANOVA_v1.csv",sep=",",col.names = NA)

set.seed(100)
perm_v2 <- adonis2(features_v2 ~ Sectio+sex+X2M_Tiere+Ern_M_vege+X2M_moisture+Straße,data=metadata_v2, permutations = 999, method = "bray", by = "margin")
write.table(perm_v2,"3.4_PERMANOVA_v2.csv",sep=",",col.names = NA)

set.seed(100)
perm_v3 <- adonis2(features_v3 ~ Sectio+sex+X4M_Tiere+Ern_M_vege+X4M_moisture+Straße+X4M_MM,data=metadata_v3, permutations = 999, method = "bray", by = "margin")
write.table(perm_v3,"3.4_PERMANOVA_v3.csv",sep=",",col.names = NA)


# ------------- MetaDeconfoundR -------------------

columns_of_interest_v1 <- c("Sectio","GA_W","age_M","KG_M","ATB_SS","PCM_SS","VBS","Toko","MB_Kontakt","Ern_M_vege","BG_M_A","BG_M_B","Geschw","Immune_Disorder_M","Rauch_M_Ex","Tiere","moisture","FN","PreHA","sex")
columns_of_interest_v2 <- c("Sectio","ATB_SS","PCM_SS","VBS","Toko","Ern_M_vege","BG_M_A","BG_M_B","Geschw","Immune_Disorder_M","X2M_respiratory","X2M_gastrointestinal","X2M_skin","X2M_vacc","Rauch_M_Ex","X2M_Tiere","X2M_moisture","X2M_FN","X2M_PreHA","X2M_Fieber","X2M_VitD","X2M_ATB_1","X2M_Hosp","X2M_Rauch_Haus","sex")
columns_of_interest_v3 <- c("Sectio","ATB_SS","PCM_SS","VBS","Toko","Ern_M_vege","BG_M_A","BG_M_B","Geschw","Immune_Disorder_M","X4M_respiratory","X4M_gastrointestinal","X4M_skin","X4M_vacc","Rauch_M_Ex","X4M_Tiere","X4M_moisture","X4M_FN","X4M_PreHA","X4M_Fieber","X4M_VitD","X4M_ATB_1","X4M_Hosp","X4M_MM","X4M_Beikost","X4M_Rauch_Haus","sex")

set.seed(100)
glom_v1_rel <- tax_glom(ps_v1_rel,taxrank = "Genus")
metadata_v1 <- as(sample_data(glom_v1_rel),"data.frame")[,columns_of_interest_v1]
sample_ids_v1 <- rownames(metadata_v1)
metadata_v1$Sectio <- as.numeric(metadata_v1$Sectio)-1
metadata_v1$sex <- as.numeric(metadata_v1$sex)-1
metadata_v1 <- apply(metadata_v1,2,function(x) as.numeric(sub(",",".",x)))
metadata_v1 <- as.data.frame(metadata_v1)
rownames(metadata_v1) <- sample_ids_v1
features_v1 <- as.data.frame(otu_table(glom_v1_rel))

metad_v1 <- MetaDeconfound(featureMat = features_v1,metaMat = metadata_v1, nnodes=6)
raw_p_v1 <- metad_v1[1]$Ps
corr_p_v1 <- metad_v1[2]$Qs
effect_size_v1 <- metad_v1[3]$Ds
status_v1 <- metad_v1[4]$status
metad_plot_v1 <- BuildHeatmap(metad_v1,featureNames = cbind(rownames(tax_table(glom_v1_rel)),tax_table(glom_v1_rel)[,"Genus"]),q_cutoff=0.05) +
  labs(title="",subtitle="",y="",x="",fill="Effect Size") + 
  theme(axis.text.y  = element_text(face="italic"))
levels(metad_plot_v1$data$metaVariableNames)[c(1,2,5)] <- c("Mother veg","Mother ex-smoker","HA babymilk")
ggsave("3.4_MetaD_Plot_v1.png",metad_plot_v1,height=5)

set.seed(100)
glom_v2_rel <- tax_glom(ps_v2_rel,taxrank = "Genus")
metadata_v2 <- as(sample_data(glom_v2_rel),"data.frame")[,columns_of_interest_v2]
sample_ids_v2 <- rownames(metadata_v2)
metadata_v2$Sectio <- as.numeric(metadata_v2$Sectio)-1
metadata_v2$sex <- as.numeric(metadata_v2$sex)-1
metadata_v2 <- apply(metadata_v2,2,function(x) as.numeric(sub(",",".",x)))
metadata_v2 <- as.data.frame(metadata_v2)
rownames(metadata_v2) <- sample_ids_v2
features_v2 <- as.data.frame(otu_table(glom_v2_rel))

metad_v2 <- MetaDeconfound(featureMat = features_v2,metaMat = metadata_v2, nnodes=6)
raw_p_v2 <- metad_v2[1]$Ps
corr_p_v2 <- metad_v2[2]$Qs
effect_size_v2 <- metad_v2[3]$Ds
status_v2 <- metad_v2[4]$status
metad_plot_v2 <- BuildHeatmap(metad_v2,featureNames = cbind(rownames(tax_table(glom_v2_rel)),tax_table(glom_v2_rel)[,"Genus"]),q_cutoff=0.05)+
  labs(title="",subtitle="",y="",x="",fill="Effect Size")+ 
  theme(axis.text.y  = element_text(face="italic"))
levels(metad_plot_v2$data$metaVariableNames) <- c("Mother ex-smoker","HA babymilk")
ggsave("3.4_MetaD_Plot_v2.png",metad_plot_v2,height=5)

set.seed(100)
glom_v3_rel <- tax_glom(ps_v3_rel,taxrank = "Genus")
metadata_v3 <- as(sample_data(glom_v3_rel),"data.frame")[,columns_of_interest_v3]
sample_ids_v3 <- rownames(metadata_v3)
metadata_v3$Sectio <- as.numeric(metadata_v3$Sectio)-1
metadata_v3$sex <- as.numeric(metadata_v3$sex)-1
metadata_v3 <- apply(metadata_v3,2,function(x) as.numeric(sub(",",".",x)))
metadata_v3 <- as.data.frame(metadata_v3)
rownames(metadata_v3) <- sample_ids_v3
features_v3 <- as.data.frame(otu_table(glom_v3_rel))

metad_v3 <- MetaDeconfound(featureMat = features_v3,metaMat = metadata_v3, nnodes=6)
raw_p_v3 <- metad_v3[1]$Ps
corr_p_v3 <- metad_v3[2]$Qs
effect_size_v3 <- metad_v3[3]$Ds
status_v3 <- metad_v3[4]$status
metad_plot_v3 <- BuildHeatmap(metad_v3,featureNames = cbind(rownames(tax_table(glom_v3_rel)),tax_table(glom_v3_rel)[,"Genus"]),q_cutoff=0.05)+
  labs(title="",subtitle="",y="",x="",fill="Effect Size")+ 
  theme(axis.text.y  = element_text(face="italic"))
levels(metad_plot_v3$data$metaVariableNames)[c(1:4,6:8)] <- c("Mother's Milk","Fever","Solid food","Siblings","Formula milk","Mother ex-smoker","HA babymilk")
ggsave("3.4_MetaD_Plot_v3.png",metad_plot_v3,height=5)

# ------------- SECTIO -----------------

# Readcounts
sig <- compare_means(DNA_ng_per_mg_Stool ~ Sectio, data = as(sample_data(ps_emma),"data.frame"),method="wilcox.test",p.adjust.method = "BH",group.by = "visit")%>% 
  mutate(fdr = format(p.adj, scientific=FALSE))%>% 
  mutate(comparison=paste(group2, group1))%>% 
  mutate(fdr.signif= case_when(p.adj <= 0.0001~"****", p.adj <= 0.001~"***", p.adj <= 0.01~"**", p.adj <= 0.049 ~ "*", p.adj>=0.05 ~ "ns"))

ggplot(sample_data(ps_emma),aes(x=Sectio,y=DNA_ng_per_mg_Stool)) + geom_boxplot(aes(fill=Sectio), show.legend = FALSE) + facet_grid(~visit) +
  scale_y_continuous(trans='log10') +
  stat_pvalue_manual(sig, label = "fdr.signif",y.position =5.4, position = "identity")+
  labs(x="Birth Mode",fill="Birth Mode",y="ng DNA / mg Stool")

ggsave("3.4_Sectio_DNA_conc.png",height=4)

# Alpha Diversity
rich_emma <- estimate_richness(ps_emma, measures = c("Shannon"))
rich_emma$visit <- sample_data(ps_emma)[rownames(rich_emma),"visit"]$visit
rich_emma$Sectio <- sample_data(ps_emma)[rownames(rich_emma),"Sectio"]$Sectio

sig <- compare_means(Shannon ~ Sectio, data = rich_emma,method="wilcox.test",p.adjust.method = "BH",group.by = "visit")%>% 
  mutate(fdr = format(p.adj, scientific=FALSE))%>% 
  mutate(comparison=paste(group2, group1))%>% 
  mutate(fdr.signif= case_when(p.adj <= 0.0001~"****", p.adj <= 0.001~"***", p.adj <= 0.01~"**", p.adj <= 0.049 ~ "*", p.adj>=0.05 ~ "ns"))

ggplot(rich_emma,aes(x=Sectio,y=Shannon)) + geom_boxplot(aes(fill=Sectio),show.legend = F) + geom_point() + facet_wrap(~visit) +
  stat_pvalue_manual(sig, label = "fdr.signif",y.position = 4, position = "identity")+
  labs(y="Alpha Diversity - Shannon Index", fill="Birth mode", x="Birth mode")

ggsave("3.4_Sectio_alpha_diversity.png",height=5)

# Beta Diversity
ps_rel.ord <- ordinate(ps_rel,"MDS","bray")
ps_v1_rel.ord <- ordinate(ps_v1_rel, "MDS", "bray") 
ps_v2_rel.ord <- ordinate(ps_v2_rel, "MDS", "bray") 
ps_v3_rel.ord <- ordinate(ps_v3_rel, "MDS", "bray") 
ord_plot_v1<- plot_ordination(ps_v1_rel, ps_v1_rel.ord, color="Sectio") + stat_ellipse(type = "t") + labs(color = "Birth mode",title="Visit 1")
ord_plot_v2<-plot_ordination(ps_v2_rel, ps_v2_rel.ord, color="Sectio") + stat_ellipse(type = "t") + labs(color = "Birth mode",title="Visit 2")
ord_plot_v3<-plot_ordination(ps_v3_rel, ps_v3_rel.ord, color="Sectio") + stat_ellipse(type = "t") + labs(color = "Birth mode",title="Visit 3")

grid.arrange(ord_plot_v1,ord_plot_v2,ord_plot_v3, ncol=1)
grid <- arrangeGrob(ord_plot_v1,ord_plot_v2,ord_plot_v3, ncol=1) 
ggsave("3.4_Sectio_beta_diversity.png", grid, height=10)

#---------- Mother milk -------------

# Overview

v1_ids <- get_variable(ps_visit1,"EMMA_ID")
v2_ids <- get_variable(ps_visit2,"EMMA_ID")
v3_ids <- get_variable(ps_visit3,"EMMA_ID")

mm_v1 <- ggplot(meta_emma_bl[meta_emma_bl$ID %in% v1_ids,],aes(x=MM)) + geom_bar(aes(fill=MM), show.legend = FALSE) + labs(x="Mother's Milk",title="Birth")  + coord_cartesian(ylim=c(0,90))
mm_v1$data$MM <- factor(mm_v1$data$MM)
levels(mm_v1$data$MM) <- c("No","Yes")

mm_v2 <- ggplot(meta_emma_2M[meta_emma_2M$ID %in% v2_ids,],aes(x=X2M_MM)) + geom_bar(aes(fill=X2M_MM), show.legend = FALSE) + labs(x="Mother's Milk",title = "2 Months") + coord_cartesian(ylim=c(0,90))
mm_v2$data$X2M_MM <- factor(mm_v2$data$X2M_MM)
levels(mm_v2$data$X2M_MM) <- c("Yes","No")
mm_v2$data$X2M_MM <- factor(mm_v2$data$X2M_MM,levels=c("No","Yes"))
mm_v2 <- mm_v2 + scale_x_discrete(drop = FALSE)+ scale_fill_manual(values=c("#00BFC4"))

mm_v3 <- ggplot(meta_emma_4M[meta_emma_4M$ID %in% v3_ids,],aes(x=X4M_MM)) + geom_bar(aes(fill=X4M_MM), show.legend = FALSE) + labs(x="Mother's Milk",title="4 Months") + coord_cartesian(ylim=c(0,90))
mm_v3$data$X4M_MM <- factor(mm_v3$data$X4M_MM) 
levels(mm_v3$data$X4M_MM) <- c("No","Yes")

grid.arrange(mm_v1,mm_v2,mm_v3, ncol=3)
grid <- arrangeGrob(mm_v1,mm_v2,mm_v3, ncol=3) 
ggsave("3.4_MM_overview.png", grid, height=3)

# Alpha diversity

# visit 3
ps_emma_4MM <- prune_samples(!is.na(sample_data(ps_emma)$X4M_MM),ps_emma)
sample_data(ps_emma_4MM)$X4M_MM <- factor(sample_data(ps_emma_4MM)$X4M_MM)
levels(sample_data(ps_emma_4MM)$X4M_MM) <- c("No Mother's Milk","Mother's Milk")

ps_visit3_4MM <- prune_samples(!is.na(sample_data(ps_visit3)$X4M_MM),ps_visit3)
sample_data(ps_visit3_4MM)$X4M_MM <- factor(sample_data(ps_visit3_4MM)$X4M_MM)
levels(sample_data(ps_visit3_4MM)$X4M_MM) <- c("No Mother's Milk","Mother's Milk")

rich_visit3 <- estimate_richness(ps_visit3_4MM, measures = "Shannon")
rich_visit3$X4M_MM <- sample_data(ps_visit3_4MM)[rownames(rich_visit3),"X4M_MM"]$X4M_MM

sig <- compare_means(Shannon ~ X4M_MM, data = rich_visit3,method="wilcox.test",p.adjust.method = "BH")%>% 
  mutate(fdr = format(p.adj, scientific=FALSE))%>% 
  mutate(comparison=paste(group2, group1))%>% 
  mutate(fdr.signif= case_when(p.adj <= 0.0001~"****", p.adj <= 0.001~"***", p.adj <= 0.01~"**", p.adj <= 0.049 ~ "*", p.adj>=0.05 ~ "ns"))

ggplot(rich_visit3,aes(x=X4M_MM,y=Shannon)) + geom_boxplot(aes(fill=X4M_MM),show.legend = F) + geom_point()  +
  stat_pvalue_manual(sig, label = "fdr.signif",y.position = 4, position = "identity")+
  labs(y="Alpha Diversity - Shannon Index", fill="Mother's milk", x="Mother's milk")

ggsave("3.4_MM_alpha_diversity_v3_Shannon.png",height=4)

# Development Shannon
rich_emma_4MM <- estimate_richness(ps_emma_4MM, measures = "Shannon")
rich_emma_4MM$X4M_MM <- sample_data(ps_emma_4MM)[rownames(rich_emma_4MM),"X4M_MM"]$X4M_MM
rich_emma_4MM$visit <- sample_data(ps_emma_4MM)[rownames(rich_emma_4MM),"visit"]$visit

sig <- compare_means(Shannon ~ visit, data = rich_emma_4MM,method="wilcox.test",p.adjust.method = "BH",group.by = "X4M_MM")%>% 
  mutate(fdr = format(p.adj, scientific=FALSE))%>% 
  mutate(comparison=paste(group2, group1))%>% 
  mutate(fdr.signif= case_when(p.adj <= 0.0001~"****", p.adj <= 0.001~"***", p.adj <= 0.01~"**", p.adj <= 0.049 ~ "*", p.adj>=0.05 ~ "ns"))

ggplot(rich_emma_4MM,aes(x=visit,y=Shannon)) + geom_boxplot(aes(fill=X4M_MM),show.legend = F) + geom_point()  + facet_wrap(~X4M_MM)+
  stat_pvalue_manual(sig, label = "fdr.signif",y.position = c(4,4.3,4.6), position = "identity")+
  labs(y="Alpha Diversity - Shannon Index", fill="Mother's milk", x="")

ggsave("3.4_MM_Shannon_development.png",height=4)

# Development Observed Species
obs_emma_4MM <- estimate_richness(ps_emma_4MM, measures = "Observed")
obs_emma_4MM$X4M_MM <- sample_data(ps_emma_4MM)[rownames(rich_emma_4MM),"X4M_MM"]$X4M_MM
obs_emma_4MM$visit <- sample_data(ps_emma_4MM)[rownames(rich_emma_4MM),"visit"]$visit

sig <- compare_means(Observed ~ visit, data = obs_emma_4MM,method="wilcox.test",p.adjust.method = "BH",group.by = "X4M_MM")%>% 
  mutate(fdr = format(p.adj, scientific=FALSE))%>% 
  mutate(comparison=paste(group2, group1))%>% 
  mutate(fdr.signif= case_when(p.adj <= 0.0001~"****", p.adj <= 0.001~"***", p.adj <= 0.01~"**", p.adj <= 0.049 ~ "*", p.adj>=0.05 ~ "ns"))

ggplot(obs_emma_4MM,aes(x=visit,y=Observed)) + geom_boxplot(aes(fill=X4M_MM),show.legend = F) + geom_point()  + facet_wrap(~X4M_MM)+
  stat_pvalue_manual(sig, label = "fdr.signif",y.position = c(96,103,110), position = "identity")+
  labs(y="Observed Species", fill="Mother's milk", x="")

ggsave("3.4_MM_Observed_development.png",height=4)

# Beta Diversity
ps_visit3_4MM_rel <- transform_sample_counts(ps_visit3_4MM, function(x) x/sum(x))
ps_visit3_4MM_rel.ord <- ordinate(ps_visit3_4MM_rel, "MDS", "bray") 
plot_ordination(ps_visit3_4MM_rel, ps_visit3_4MM_rel.ord, color="X4M_MM") + stat_ellipse(type = "t") + labs(color = "Mother's milk")
ggsave("3.4_MM_beta_diversity_v3.png",height = 4)
