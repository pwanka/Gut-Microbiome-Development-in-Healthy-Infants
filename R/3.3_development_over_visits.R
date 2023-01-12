library(phyloseq)
library(ggplot2)
library("RColorBrewer")
library(ggpubr)
library(dplyr)
library(vegan)

emma_rds_file <- "ps_emma.rds"
ps_emma <- readRDS(emma_rds_file)

ps_rel <- transform_sample_counts(ps_emma,function(x) x/sum(x))

ps_visit1 <- prune_samples(get_variable(ps_emma, "visit") == "Visit 1" , ps_emma)
ps_visit2 <- prune_samples(get_variable(ps_emma, "visit") == "Visit 2" , ps_emma)
ps_visit3 <- prune_samples(get_variable(ps_emma, "visit") == "Visit 3" , ps_emma)

# Make nice Barplot
nice_rel_barplot <- function(ps_object,taxrank="Phylum",minAbundance=0.02){
  
  ps <- tax_glom(ps_object,taxrank = taxrank) 
  ps <- transform_sample_counts(ps, function(x) x / sum(x))
  
  overmin<- colSums(otu_table(ps) > minAbundance)
  for (name in names(overmin[overmin == 0])){
    tax_table(ps)[name,taxrank] <- paste("Abundance < ",minAbundance*100,"% in all Samples",sep="")
  }

  plot <- plot_bar(ps, fill = taxrank) + geom_bar(aes_string(color = taxrank, fill=taxrank), stat="identity", position="stack") + 
    theme(strip.text.x = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank())
  
  sorted <- plot$data[!(duplicated(plot$data$Sample)),]
  sorted <- sorted[with(sorted,order(sorted$Class,sorted$Abundance)),]
  sampleorder <- sorted$Sample
  
  plot$data$Sample <- factor(plot$data$Sample, levels = sampleorder)
  return(plot)
}

# make Palette to have same colors for same taxa across different plots
getPalette = colorRampPalette(brewer.pal(8, "Set1"))
phylaList = unique(tax_table(ps_emma)[,"Class"])
phylaPalette = getPalette(length(phylaList))
names(phylaPalette) = phylaList

# taxonomy barplots
nice_rel_barplot(ps_visit1,taxrank = "Class") + 
  scale_fill_manual(values=phylaPalette,limits=force) +
  scale_color_manual(values=phylaPalette,limits=force) + 
  labs(y="Relative Abundance", x = "Visit 1 Samples")
ggsave("3.3_Barplot_Class_v1.png")

nice_rel_barplot(ps_visit2,taxrank="Class") + 
  scale_fill_manual(values=phylaPalette,limits=force) +
  scale_color_manual(values=phylaPalette,limits=force) + 
  labs(y="Relative Abundance", x = "Visit 2 Samples")
ggsave("3.3_Barplot_Class_v2.png")

nice_rel_barplot(ps_visit3,taxrank="Class") + 
  scale_fill_manual(values=phylaPalette,limits=force) +
  scale_color_manual(values=phylaPalette,limits=force) + 
  labs(y="Relative Abundance", x = "Visit 3 Samples") 
ggsave("3.3_Barplot_Class_v3.png")

# Top 5 Genera
genera_all <- tax_glom(ps_rel,taxrank = "Genus")

genera_v1 <- prune_samples(get_variable(genera_all, "visit") == "Visit 1" , genera_all)
top5_median_v1 <- sort(apply(otu_table(genera_v1),2,median),decreasing = T)[1:5]
top5_mean_v1 <- colMeans(otu_table(genera_v1)[,names(top5_median_v1)])
top5_prevalence_v1 <- colSums(otu_table(genera_v1)[,names(top5_median_v1)]>0)/nrow(otu_table(genera_v1)[,names(top5_median_v1)])
top5_v1 <- cbind(top5_median_v1,top5_mean_v1,top5_prevalence_v1)
rownames(top5_v1) <- tax_table(genera_v1)[rownames(top5_v1),"Genus"]
colnames(top5_v1) <- c("Median","Mean","Prevalence")
write.table(top5_v1,"3.3_Top5_Genera_v1.csv",sep=",",col.names = NA)

genera_v2 <- prune_samples(get_variable(genera_all, "visit") == "Visit 2" , genera_all)
top5_median_v2 <- sort(apply(otu_table(genera_v2),2,median),decreasing = T)[1:5]
top5_mean_v2 <- colMeans(otu_table(genera_v2)[,names(top5_median_v2)])
top5_prevalence_v2 <- colSums(otu_table(genera_v2)[,names(top5_median_v2)]>0)/nrow(otu_table(genera_v2)[,names(top5_median_v2)])
top5_v2 <- cbind(top5_median_v2,top5_mean_v2,top5_prevalence_v2)
rownames(top5_v2) <- tax_table(genera_v2)[rownames(top5_v2),"Genus"]
colnames(top5_v2) <- c("Median","Mean","Prevalence")
write.table(top5_v2,"3.3_Top5_Genera_v2.csv",sep=",",col.names = NA)

genera_v3 <- prune_samples(get_variable(genera_all, "visit") == "Visit 3" , genera_all)
top5_median_v3 <- sort(apply(otu_table(genera_v3),2,median),decreasing = T)[1:5]
top5_mean_v3 <- colMeans(otu_table(genera_v3)[,names(top5_median_v3)])
top5_prevalence_v3 <- colSums(otu_table(genera_v3)[,names(top5_median_v3)]>0)/nrow(otu_table(genera_v3)[,names(top5_median_v3)])
top5_v3 <- cbind(top5_median_v3,top5_mean_v3,top5_prevalence_v3)
rownames(top5_v3) <- tax_table(genera_v3)[rownames(top5_v3),"Genus"]
colnames(top5_v3) <- c("Median","Mean","Prevalence")
write.table(top5_v3,"3.3_Top5_Genera_v3.csv",sep=",",col.names = NA)

genera_to_compare <- unique(c(names(top5_median_v1),names(top5_median_v2),names(top5_median_v3)))
comparison_data <- psmelt(ps_rel)
comparison_data <- comparison_data[comparison_data$OTU %in% genera_to_compare,]
comparisons <- list(c("Visit 1", "Visit 2"), c("Visit 1", "Visit 3"),c("Visit 2", "Visit 3"))

sig <- compare_means(Abundance ~ visit, data = comparison_data,method="wilcox.test",p.adjust.method = "BH",group.by = "Genus")%>% 
  mutate(fdr = format(p.adj, scientific=FALSE))%>% 
  mutate(comparison=paste(group2, group1))%>% 
  mutate(fdr.signif= case_when(p.adj <= 0.0001~"****", p.adj <= 0.001~"***", p.adj <= 0.01~"**", p.adj <= 0.049 ~ "*", p.adj>=0.05 ~ "ns"))

ggplot(comparison_data,aes(x=visit,y=Abundance)) + geom_boxplot(aes(fill=visit)) + facet_wrap(~Genus,nrow=2) +
  labs(y="Relative Abundance",x="Visit",fill="Visit") +
  stat_pvalue_manual(sig, label = "fdr.signif",y.position = c(1.1,1.25,1.4), position = "identity") +                                                             # Change font size
  theme(strip.text.x = element_text(size = 7)) 

ggsave("3.3_Top_Genera_over_visits.png",height=8)

rich_emma <- estimate_richness(ps_emma, measures = c("Shannon"))
rich_emma$visit <- sample_data(ps_emma)[rownames(rich_emma),"visit"]$visit

sig <- compare_means(Shannon ~ visit, data = rich_emma,method="wilcox.test",p.adjust.method = "BH")%>% 
  mutate(fdr = format(p.adj, scientific=FALSE))%>% 
  mutate(comparison=paste(group2, group1))%>% 
  mutate(fdr.signif= case_when(p.adj <= 0.0001~"****", p.adj <= 0.001~"***", p.adj <= 0.01~"**", p.adj <= 0.049 ~ "*", p.adj>=0.05 ~ "ns"))

ggplot(rich_emma,aes(x=visit,y=Shannon)) + geom_boxplot(aes(fill=visit)) + geom_point() +
  stat_pvalue_manual(sig, label = "fdr.signif",y.position = c(4,4.3,4.6), position = "identity")+
  labs(y="Alpha Diversity - Shannon Index", fill="Visit", x="Visit")
ggsave("3.3_Alpha_Diversity_Visits.png",height=5)

# Beta Diversity
ps_rel.ord <- ordinate(ps_rel, "MDS", "bray") 
plot_ordination(ps_rel, ps_rel.ord, color="visit") + stat_ellipse(type = "t") + labs(color="Visit")
ggsave("3.3_Beta_diversity_visits.png",height=6)

