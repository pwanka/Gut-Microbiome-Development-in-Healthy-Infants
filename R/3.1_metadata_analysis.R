library(phyloseq)
library(ggplot2)
library(gridExtra)

all_rds_file <- "ps_all.rds"
emma_rds_file <- "ps_emma_all_samples.rds"
ps_all <- readRDS(all_rds_file)
ps_emma <- readRDS(emma_rds_file)
metadata <- sample_data(ps_emma)

emma_metafile_baseline <- "Data/220214_EMMA Datenbank-BL.csv"
meta_emma_bl <- read.csv(emma_metafile_baseline, sep = "\t", skip = 2)[-1,]

# Barplot sex + birthmode
meta_emma_bl$sex <- factor(meta_emma_bl$sex)
levels(meta_emma_bl$sex)[1] = "Female"
levels(meta_emma_bl$sex)[2] = "Male"
plot_sex <- ggplot(meta_emma_bl, aes(sex,fill=sex)) + geom_bar(show.legend = F) +
  labs(title="A",x="Sex", fill ="Sex",y="")

meta_emma_bl$Sectio <- factor(meta_emma_bl$Geb_Modus %in% c(1,2,3))
levels(meta_emma_bl$Sectio)[1] = "Vaginal"
levels(meta_emma_bl$Sectio)[2] = "Sectio"
plot_sectio <- ggplot(meta_emma_bl, aes(Sectio,fill=Sectio)) + geom_bar(show.legend = F) +
  labs(title="B", x="Birth Mode", fill = "Birth Mode",y="")

grid.arrange(plot_sex, plot_sectio, ncol=2)
grid <- arrangeGrob(plot_sex, plot_sectio, ncol=2) 
ggsave("3.1_sex_birthmode.png", grid,height=4)

# Readcounts per visit 

ggplot(metadata,aes(x=visit,y=total_reads,fill=visit)) + geom_boxplot() +
  scale_y_continuous(trans='log10') +
  labs(x="Visit", fill = "Visit",y="Readcounts") +
  geom_hline(yintercept=1000,linetype="dashed",color="red") 
ggsave("3.1_Readcounts_per_visit.png",height=6)

# DNA content: viel na :( --> fill in DNA NA values from replicates of same sample in different run
dnana <- metadata$X.SampleID[is.na(metadata$DNA_ng_per_mg_Stool)] #get samples with na in dna column
reps <- sample_data(ps_all)[substr(sample_data(ps_all)$X.SampleID,nchar(sample_data(ps_all)$X.SampleID)-6,nchar(sample_data(ps_all)$X.SampleID)) %in% substr(dnana,nchar(dnana)-6,nchar(dnana)) & !is.na(sample_data(ps_all)$DNA_ng_per_mg_Stool),c("X.SampleID","DNA_ng_per_mg_Stool")]
rownames(reps)<- substr(reps$X.SampleID,nchar(reps$X.SampleID)-6,nchar(reps$X.SampleID))

metadata_new <- as(metadata,"data.frame") %>% 
  mutate(DNA_ng_per_mg_Stool = ifelse(is.na(DNA_ng_per_mg_Stool),as(reps,"data.frame")[substr(X.SampleID,nchar(X.SampleID)-6,nchar(X.SampleID)),"DNA_ng_per_mg_Stool"],DNA_ng_per_mg_Stool))

# DNA content in samples Plot + Summary Table
plot_dna_per_visit <- ggplot(subset(metadata_new,DNA_ng_per_mg_Stool < 1000),aes(x=visit,y=DNA_ng_per_mg_Stool,fill=visit)) + geom_boxplot() + 
  labs(x="Visit", fill = "Visit", y="ng DNA / mg Stool",title="A") 
# 2 datapoints from visit 3 removed for legibility

sv1 <- summary(metadata_new$DNA_ng_per_mg_Stool[metadata_new$visit=="Visit 1"])
sv2 <- summary(metadata_new$DNA_ng_per_mg_Stool[metadata_new$visit=="Visit 2"])
sv3 <- summary(metadata_new$DNA_ng_per_mg_Stool[metadata_new$visit=="Visit 3"])

dna_content_table <- round(rbind(sv1,sv2,sv3)[,1:6],2)
rownames(dna_content_table) <- paste("Visit",c(1:3));dna_content_table
write.table(dna_content_table,"DNA_content_Table.csv",sep=",")


# Samples <1ngDNA/mgStool
metadata_new$lowDNA <- factor(metadata_new$DNA_ng_per_mg_Stool < 1)
levels(metadata_new$lowDNA)[1] = "≥ 1ng DNA / mg Stool"
levels(metadata_new$lowDNA)[2] = "< 1ng DNA / mg Stool"
plot_low_dna <- ggplot(metadata_new,aes(x=visit,fill=lowDNA)) + geom_bar() +
labs(x="Visit", fill = "DNA content", y="",title = "B") 

grid.arrange(plot_dna_per_visit, plot_low_dna, ncol=1)
grid <- arrangeGrob(plot_dna_per_visit, plot_low_dna, ncol=1) 
ggsave("3.1_DNA_content_per_visit.png", grid, height=6)

# Readcount vs DNA concentration
ggplot(subset(metadata_new,DNA_ng_per_mg_Stool < 100),aes(x=DNA_ng_per_mg_Stool,y=total_reads)) + geom_point()+ 
  labs(x="ng DNA / mg Stool",y="Readcount") 
ggsave("3.1_readcount_vs_dna.png",height=4)
cor(subset(metadata_new,DNA_ng_per_mg_Stool < 100)$DNA_ng_per_mg_Stool,subset(metadata_new,DNA_ng_per_mg_Stool < 100)$total_reads,use = "complete.obs")

# Readcount vs diversity (by visit)
div <- estimate_richness(ps_emma,measures = "Shannon")
metadata$Shannon <- div$Shannon[match(rownames(metadata),rownames(div))] 

ggplot(metadata,aes(x=total_reads,y=Shannon)) + geom_point(aes(color=visit),show.legend = F) + facet_grid(~visit) + 
  labs(x="Readcount", y="Alpha Diversity - Shannon Index") 
ggsave("3.1_readcount_vs_diversity.png",height=4)
