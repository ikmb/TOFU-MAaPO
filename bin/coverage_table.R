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

contigsTPMtable <- data.table(contigName=depths$contigName, TPMpercontig, contigLen=depths$contigLen)

## Solution with mean avg per contig
#output <- merge(binstaxcontigs, contigsTPMtable, by.x="contig", by.y="contigName") %>% group_by(binnew) %>% summarise(mean_contigs_TPM_per_bin=(mean(TPMpercontig)),taxa=unique(classification)) %>% rename(bin=binnew)

## Solution with weighted by contig_length
output <- merge(binstaxcontigs, contigsTPMtable, by.x="contig", by.y="contigName") %>% 
  group_by(binnew) %>% 
  summarise(weighted_contigs_TPM_per_bin=sum((TPMpercontig*contigLen)/sum(contigLen)),
            taxa=unique(classification)) %>% 
  summarise(bin=binnew,
            taxa=binstotax$classification,
            scaled_weighted_contigs_TPM_per_bin=weighted_contigs_TPM_per_bin/sum(weighted_contigs_TPM_per_bin)*10^6) 

fwrite(output, file=output_file_name,sep = ",")
