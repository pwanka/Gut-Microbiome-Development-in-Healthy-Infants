#!/usr/bin/env Rscript

rm(list=ls())


library("optparse")
option_list = list(
  make_option(c("-p", "--inpath"), type="character", default=NULL, 
              help="path were metadata and sequence files are stored", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL, 
              help="path to write files to", metavar="character"),
  make_option(c("-r", "--runnum"), type="character", default=NULL, 
              help="MiSeq sequencing run number", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$inpath) | is.null(opt$outpath) | is.null(opt$runnum)){
  print_help(opt_parser)
  stop("please supply the necessary arguments", call.=FALSE)
}


start_time = Sys.time()

# config
# put in paths to the raw sequence-files and the metadata 
path <- opt$inpath
setwd(path)
out_path <- opt$outpath
run_num <- opt$runnum

metadata_file = file.path(path,paste(run_num, "Sample_Metadata.csv", sep="_"))
# this one is just used for the names of the outputfiles
#--------
library(dada2)
packageVersion("dada2")
library("phyloseq")
packageVersion("phyloseq")
run_num <- as.character(run_num)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
#plotQualityProfile(fnFs)
#plotQualityProfile(fnRs)

filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# apply filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, truncLen = c(240,240), trimLeft = c(20,20), rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# filter out files were quality trimming removed all reads
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]
out
# subset the out object to only contain files that did not end up with 0 reads after filtering
out2 <- as.data.frame(out)
out_filter <- out2$reads.out != 0
out2 <- out2[out_filter,]
out2 <- as.matrix(out2)
out <- out2

# learn error model
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# adapt sample names to filtered reads (when some might have been removed due to filtering)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.names

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
table(seqtab)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#cutting a band in silico
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(249,254)]
seqtab <- seqtab2

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# percent of sequences that are not chimeric
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track
track_file_name <- paste(run_num,"Sample_DADA_statistics.csv", sep="_")
track_file <- paste(out_path,track_file_name,sep="/")
write.table(track, file = track_file, quote = FALSE, sep = "\t", col.names=NA)

session_name <- paste(out_path,paste(run_num,"R_session_DADA2_before_taxonomy_assignment.Rdata", sep="_"),sep="/")
save.image(session_name)

print("assigning taxonomy")
# SILVA Assignment V.132

taxa <- assignTaxonomy(seqtab.nochim, "Data/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa <- addSpecies(taxa, "Data/silva_species_assignment_v138.1.fa.gz")


# RDP Assignemnt

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print,30)

session_name <- paste(out_path,paste(run_num,"R_session_DADA2.Rdata", sep="_"),sep="/")
save.image(session_name)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))

otu_output_file <- paste(out_path,paste(run_num,"ASV_table.csv", sep="_"),sep="/")
tax_output_file <- paste(out_path,paste(run_num,"ASV_taxonomy.csv", sep="_"),sep="/")
otus <- otu_table(ps)
tax <- tax_table(ps)
write.table(otus, file = otu_output_file, quote = FALSE, sep = "\t", col.names=NA)
write.table(tax, file = tax_output_file, quote = FALSE, sep = "\t", col.names=NA)

session_name <- paste(out_path,paste(run_num,"R_session_DADA2.Rdata", sep="_"),sep="/")
save.image(session_name)
ps_output_file <- paste(out_path,paste(run_num,"DADA2_phyloseqobject.rds", sep="_"),sep="/")
saveRDS(ps , ps_output_file)

#------------------------------------- DECONTAM ------------------------------------------------------------
library("tidyverse")
library("decontam")

metadata <- read.csv(metadata_file, fill = TRUE, sep = "\t")
metadata <- metadata[!metadata$X.SampleID == "",]
rownames(metadata) <- metadata$X.SampleID
sample_data <- sample_data(metadata)
ps <- merge_phyloseq(ps, sample_data)
str(sample_data(ps)$dna_quant_ng_ul)
ps.clean <- subset_taxa(ps, !Kingdom %in% c("Archaea","Eukaryota") & !is.na(Phylum) & !Phylum %in% c("","uncharacterized"))
ps <- ps.clean
#--------- update the run statistics --------------
track_update <- cbind(track, rowSums(otu_table(ps.clean)))
colnames(track_update) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim","NOarchaeaEukaryota")
rownames(track_update) <- sample.names
track_update
write.table(track_update, file = track_file, quote = FALSE, sep = "\t", col.names=NA)
#--------------------------------------------------
saveRDS(ps , ps_output_file)

sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
# remove samples with 0 read counts (samples with fewer reads than defined in minreads)
minreads = 1
if (any(sample_sums(ps) < minreads)) {
  ps.pre <- prune_samples(sample_sums(ps) < minreads, ps)
  sample_names(ps.pre)
} else {cat(paste("there were no samples with less than", minreads,"read", sep=" "))}

ps <- prune_samples(sample_sums(ps) >= minreads, ps)

# set 0 and NA in the DNA quantification measurements performed by qubit 4 high sensetivity (HS) to 0.0005 ng/ul (minimal resolution) 
# this is necessary, because decontam can't calculate DNA quantifications == 0
sample_data(ps)$dna_quant_ng_ul[sample_data(ps)$dna_quant_ng_ul == 0] = 0.0005
sample_data(ps)$dna_quant_ng_ul[is.na(sample_data(ps)$dna_quant_ng_ul)] = 0.0005
contamdf.com <- isContaminant(ps,conc = "dna_quant_ng_ul", neg = "is.neg", method = "combined", threshold = c(0.05))

table(contamdf.com$contaminant)
head(which(contamdf.com$contaminant))
ps.noncontam <- prune_taxa(!contamdf.com$contaminant, ps)

otu_output_file <- paste(out_path,paste(run_num,"ASV_table_decontam.csv", sep="_"),sep="/")
tax_output_file <- paste(out_path,paste(run_num,"ASV_taxonomy_decontam.csv", sep="_"),sep="/")

otus <- otu_table(ps.noncontam)
tax <- tax_table(ps.noncontam)
ps.noncontam
write.table(otus, file = otu_output_file, quote = FALSE, sep = "\t", col.names=NA)
write.table(tax, file = tax_output_file, quote = FALSE, sep = "\t", col.names=NA)

session_name <- paste(out_path,paste(run_num,"R_session_DADA2_decontam.Rdata", sep="_"),sep="/")
save.image(session_name)
ps_output_file <- paste(out_path,paste(run_num,"DADA2_phyloseqobject_decontam.rds",sep="_"),sep="/")
saveRDS(ps , ps_output_file)

end_time = Sys.time()
start_time
end_time
end_time-start_time
