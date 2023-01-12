library(phyloseq)
library("ggplot2")
library(dplyr)
library(gridExtra)

nks<- c("NK 04.05.2021","PCR NK 04.05.2021","NK 15.07.2021","PCR NK 15.07.2021","NK 21.07.2021","PCR NK 21.07..2021","NK 12.08.2021","NK 16.08.2021","NK 18.08.2021","NK 20.08.2021","NK PCR 12.08.21","NK PCR 16.08.21","NK PCR 18.08.21","NK PCR 20.08.21","NK 09.11.2021","NK PCR 09.11.21","NK 21.09.2021","NK 23.09.21","NK 25.09.2021","NK 04.10.2021","NK 01.10.2021","NK 21.07.2021","NK 25.10.2021","NK PCR 25.10.21","NK PCR 21.09.21","NK PCR 23.09.21","NK PCR 25.09.21","NK PCR 04.10.21","NK PCR 01.10.21","NK PCR 21.07.21","NK 27.10.2021","NK 12.11.2021","NK PCR 12.11.21","NK PCR 27.10.2021")

all_rds_file <- "ps_all.rds"
ps_all <- readRDS(all_rds_file)

emma_rds_file <- "ps_emma.rds"
ps_emma <- readRDS(emma_rds_file)

ps_nk <- subset_samples(ps_all, material == "Water")
#ps_nk <- subset_samples(ps_nk, run %in% unique(get_variable(ps_emma,"run")))
ps_nk <- subset_samples(ps_nk, Label %in% nks)
ps_nk <- prune_taxa(taxa_sums(ps_nk) > 0, ps_nk)

sample_data(ps_nk)$total_reads <- sample_sums(ps_nk)

# boxplot readcounts
plot_box <- ggplot(sample_data(ps_nk),aes(y=total_reads)) + 
  geom_boxplot(fill="orange") +
  theme(axis.text.y = element_blank(),axis.ticks = element_blank()) +
  labs(title="A",y= "Readcount") + 
  geom_hline(yintercept=1000,linetype="dashed")  +
  coord_flip()

# taxonomy 

plot_tax_nk <- plot_bar(ps_nk,fill="Class") +
  geom_bar(aes(fill=Class), stat="identity", position="stack") + 
  theme(axis.text.x = element_blank(),axis.ticks = element_blank())+ 
  labs(title= "B",y= "Readcount",x="") + 
  geom_hline(yintercept=1000,linetype="dashed") 

grid.arrange(plot_box, plot_tax_nk, ncol=1)
grid <- arrangeGrob(plot_box, plot_tax_nk, ncol=1) 
ggsave("3.2_nk_boxplot_tax.png", grid)

# compare >1000 reads with affected samples
ps_many <- prune_samples(sample_sums(ps_nk)>=1000, ps_nk)

affected_dates <- get_variable(ps_many,"extraction_date")
ps_affected <- subset_samples(ps_emma, extraction_date %in% affected_dates)
ps_affected_nk <- merge_phyloseq(ps_many,ps_affected)
sample_data(ps_affected_nk)$nk <- factor(get_variable(ps_affected_nk,"Label") %in% nks)
levels(sample_data(ps_affected_nk)$nk)[1] = ""
levels(sample_data(ps_affected_nk)$nk)[2] = "C"
sample_data(ps_affected_nk)$extraction_date <- factor(sample_data(ps_affected_nk)$extraction_date)

p <- plot_bar(ps_affected_nk,fill="Class") +
  geom_bar(aes(fill=Class), stat="identity", position="stack") + 
  labs(y= "Readcount") + 
  facet_wrap(~extraction_date, scales="free_x", nrow=1)+ 
  theme(strip.text.x = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank())+
  geom_text(aes(label=nk, y=0),vjust=2,size=2)
p$data$Class[p$data$Class == ""] <- "Undefined"
ggsave("3.2_nk_vs_affected.png",p)

# ordination plot nach extraction date mit relative abundances
levels(sample_data(ps_rel_aff)$nk)[1] = "Sample"
levels(sample_data(ps_rel_aff)$nk)[2] = "Negative Control"
levels(sample_data(ps_rel_aff)$extraction_date) <- paste("Extraction Date",1:4)
ps_rel_affected_nk.ord <- ordinate(ps_rel_aff, "MDS", "bray")
plot_ordination(ps_rel_aff, ps_rel_affected_nk.ord, color="extraction_date",shape="nk") + labs(color="Extraction Date", shape="Sample/Control")
ggsave("3.2_ordination_affected_nk.png")
