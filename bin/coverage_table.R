#!/usr/bin/env Rscript
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly=TRUE)
depths <- fread(args[1])
expected_avg_read_length <- as.integer(args[2])
contstobins <- fread(args[3])
classification <- fread(args[4])
sampleID <- args[5]
output_file_name <- args[6]

binstotax <- classification %>% select(user_genome, classification)
binstaxcontigs <- merge(contstobins, binstotax, by.x="binnew", by.y="user_genome")

depths$contigLen <- as.numeric(depths$contigLen)
depths$totalAvgDepth <- as.numeric(depths$totalAvgDepth)

expected_transcripts_percontig <- (depths$contigLen*depths$totalAvgDepth)/expected_avg_read_length

#TPM formula
#reads mapped to transcript divided by length of transcript divided by the sum of all reads mapped to transcript divided by length of transcript multiplied by one million

TPMpercontig <- (expected_transcripts_percontig/depths$contigLen)/sum((expected_transcripts_percontig/depths$contigLen))*10^6

contigsTPMtable <- data.table(contigName=depths$contigName, TPMpercontig)



output <- merge(binstaxcontigs, contigsTPMtable, by.x="contig", by.y="contigName") %>% group_by(binnew) %>% summarise(mean_contigs_TPM_per_bin=(mean(TPMpercontig)),taxa=unique(classification)) %>% rename(bin=binnew)

fwrite(output, file=output_file_name,sep = ",")
#paste0(sampleID,"_abundance_table.tbl")

#
