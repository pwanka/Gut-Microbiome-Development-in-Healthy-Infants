library(phyloseq)

asvfile <- "Data/ASV_table_decontam_combined.csv"
taxfile <- "Data/ASV_taxonomy_decontam_combined.csv"
metafile <- "Data/Metadata_combined_adjusted.csv"
emma_metafile_baseline <- "Data/220214_EMMA Datenbank-BL.csv"
emma_metafile_2M <- "Data/220214_EMMA Datenbank-2M.csv"
emma_metafile_4M <- "Data/220214_EMMA Datenbank-4M.csv"

asv <- otu_table(read.csv(asvfile, sep = "\t", row.names = 1),taxa_are_rows = FALSE)
tax <- read.csv(taxfile, sep = "\t", row.names = 1)
tax_matrix <- as(tax, 'matrix')
tax <- tax_table(tax_matrix)
metadata <- read.csv(metafile, fill = TRUE, sep = "\t")
metadata <- subset(metadata, X.SampleID != "")
metadata <- metadata[-c(3091:3097),] #remove duplicate rows
rownames(metadata) <- metadata$X.SampleID
sample_data <- sample_data(metadata)

ps_all <- phyloseq(asv, tax, sample_data)

dna <- Biostrings::DNAStringSet(taxa_names(ps_all))
names(dna) <- taxa_names(ps_all)
ps_all <- merge_phyloseq(ps_all, dna)
taxa_names(ps_all) <- paste0("ASV", seq(ntaxa(ps_all)))

# add run column
runs<- c()
for (run in substr(sample_data(ps_all)$X.SampleID,1,2))
{
  pruns <- as.character(1:20)
  if (!(run %in% pruns))
  {
    run <- substr(run,1,1)
  }
  runs <- append(runs,run)
}
sample_data(ps_all)$run <- runs

ps_all_output_file <- "ps_all.rds"
saveRDS(ps_all, ps_all_output_file)

ps_emma <- subset_samples(ps_all, project == "EMMA")
ps_emma <- subset_samples(ps_emma, exclusion != "True")
ps_emma <- prune_taxa(taxa_sums(ps_emma) > 0, ps_emma)

# rename mistakes in samplenames
sample_data(ps_emma)$Label[sample_data(ps_emma)$Label == "15012S021c002"] <- "01.2.S.021.c.002"
sample_data(ps_emma)$Label[sample_data(ps_emma)$Label == "15012V022c028"] <- "01.2.V.022.c.028"
sample_data(ps_emma)$Label[sample_data(ps_emma)$Label == "15012V033c001"] <- "01.2.V.033.c.001"
sample_data(ps_emma)$Label[sample_data(ps_emma)$Label == "01.2.V.025c028"] <- "01.2.V.015.c.028"
sample_data(ps_emma)$Label[sample_data(ps_emma)$Label == "1.1.V.057.c.090"] <- "01.1.V.057.c.090"
sample_data(ps_emma)$Label[sample_data(ps_emma)$Label == "1.2.V.018.c.028"] <- "01.2.V.018.c.028"

#remove negative controls
ps_emma <- prune_samples(substr(sample_data(ps_emma)$Label,1,2) != "NK",ps_emma)

#get EMMA ID from Sample Label
sample_data(ps_emma)$EMMA_ID <- substr(sample_data(ps_emma)$Label,1,10)

#Merge Metadata sheets to one sheet
meta_emma_bl <- read.csv(emma_metafile_baseline, sep = "\t", dec = ",", skip = 2)[-1,]
meta_emma_2M <- read.csv(emma_metafile_2M, sep = "\t",dec = ",", skip = 2)[-1,]
meta_emma_2M$ID[meta_emma_2M$ID == "01.2.S.030"] <- "01.1.S.030"
meta_emma_4M <- read.csv(emma_metafile_4M, sep = "\t",dec = ",", skip = 2)[-1,]
meta_emma_4M$ID[meta_emma_4M$ID == "01.2.S.030"] <- "01.1.S.030"
meta_emma <- merge(meta_emma_bl,meta_emma_2M, by = "ID", all.x = T)
meta_emma <- merge(meta_emma,meta_emma_4M, by = "ID", all.x = T)

colnames(meta_emma)[1] <- "EMMA_ID"

#Merge Emma-Metadata with Sequencing Metadata and add to object 
kk <- list(sample_data(ps_emma))
mergedmeta <- merge(kk,meta_emma, by = "EMMA_ID", all.x = T)
rownames(mergedmeta) <- mergedmeta$X.SampleID
sample_mergedmeta <- sample_data(mergedmeta)
sample_data(ps_emma) <- sample_mergedmeta

sample_data(ps_emma)$patient_id <- factor(sapply(strsplit(as.character(get_variable(ps_emma, "Label")), ".",fixed = TRUE), `[`, 4))
sample_data(ps_emma)$day <- factor(sapply(strsplit(as.character(get_variable(ps_emma, "Label")), ".",fixed = TRUE), `[`, 6))
sample_data(ps_emma)$Sectio <- factor(get_variable(ps_emma,"Geb_Modus") %in% c(1,2,3))
sample_data(ps_emma)$sex <- factor(sapply(strsplit(as.character(get_variable(ps_emma, "Label")), ".",fixed = TRUE), `[`, 2))
levels(sample_data(ps_emma)$sex)[1] <- "m" 
levels(sample_data(ps_emma)$sex)[2] <- "f"
sample_data(ps_emma)$total_reads <- sample_sums(ps_emma)

sample_data(ps_emma)$visit <- sample_data(ps_emma)$day
levels(sample_data(ps_emma)$visit)[1:2] <- "Visit 1"
levels(sample_data(ps_emma)$visit)[2:3] <- "Visit 2"
levels(sample_data(ps_emma)$visit)[3] <- "Visit 3"

saveRDS(ps_emma, "ps_emma_all_samples.rds")

# remove Samples with readcount < 1000 (see 3.1)
ps_emma <- subset_samples(ps_emma, total_reads >= 1000)
saveRDS(ps_emma, "ps_emma.rds")

write.table(otu_table(ps_emma), file = "Data/EMMA_ASV_table.csv", quote = FALSE, sep = "\t", col.names=NA)
write.table(tax_table(ps_emma), file = "Data/EMMA_ASV_taxonomy.csv", quote = FALSE, sep = "\t", col.names=NA)
write.table(sample_data(ps_emma), file = "Data/EMMA_Sample_Metadata.csv", quote = FALSE, sep = "\t", col.names=NA)
