# DADA2
library(dada2); packageVersion("dada2")
path <- "datos 2" 
list.files(path)

fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnFs
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
fnRs
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtFs
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
filtRs
names(filtFs) <- sample.names
names(filtFs)
names(filtRs) <- sample.names
names(filtRs)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,compress=TRUE,multithread=TRUE)
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# 0.4904637-> quimeras 

getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), 
              sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v128_train_set.fa.gz", 
                       multithread=TRUE)
taxa <- addSpecies(taxa, "silva_species_assignment_v128.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)



library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(),
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps



plot_richness(ps, measures=c("Shannon", "Simpson")) 

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, title="Bray NMDS")


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:70]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Genus")
plot_bar(ps.top20, fill="Family") 
plot_bar(ps.top20, fill="Phylum")
plot_bar(ps.top20, fill="Class")
plot_bar(ps.top20, fill="Order")
plot_bar(ps.top20, fill="Kingdom")


plot_heatmap(ps.top20)
plot_heatmap(ps.top20, "NMDS", "jaccard")

heatmap(otu_table(ps.top20))
