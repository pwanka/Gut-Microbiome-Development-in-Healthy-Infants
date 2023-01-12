library(phyloseq)
library(scales)

emma_rds_file <- "ps_emma_all_samples.rds"
ps_emma <- readRDS(emma_rds_file)

# Read filtering
filtervars <- c("input","filtered","merged","tabled","nonchim","NOarchaeaEukaryota","total_reads")
reads_filtered <- get_variable(ps_emma,filtervars)/get_variable(ps_emma,"input")
reads_filtered_mean <- sapply(colMeans(reads_filtered,na.rm=T),percent,accuracy=0.01) 
reads_filtered_iqr <- apply(apply(reads_filtered,2,quantile,probs=c(0.25,0.75),na.rm=T),2,percent,accuracy=0.01) 
reads_filtered_iqr <- paste("(",reads_filtered_iqr[1,],"-",reads_filtered_iqr[2,],")",sep="")
filtering_steps <- paste(reads_filtered_mean,reads_filtered_iqr)
names(filtering_steps) <- c("Raw reads","Quality-filtered","Merged","Correct length","Non-chimeric","Phylum-assigned Bacteria","Decontaminated")
filtering_steps <-t(as.data.frame(filtering_steps))
rownames(filtering_steps) <- "Mean (IQR) remaining % of raw reads"

write.table(filtering_steps,"2.3_filtering_steps.csv",sep=",",col.names = NA)

