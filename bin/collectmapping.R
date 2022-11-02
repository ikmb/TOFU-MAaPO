#!/usr/bin/env Rscript
##### By Malte Rühlemann & Eike Matthias Wacker 2022 #####

library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
alldepth <- args[1]
files=list.files()
samples=files[grep(".depth.txt", files)]


all_abu = sapply(samples, function(s) read.table(s, head=T, stringsAsFactors=F), simplify=F)

for(i in 1:length(all_abu)){if(i==1){all_abu_frame=all_abu[[1]]; next}; all_abu_frame=all_abu_frame %>% left_join(all_abu[[i]] %>% select(-totalAvgDepth))}

all_abu_frame$totalAvgDepth = data.frame(all_abu_frame[,grep(".bam$", colnames(all_abu_frame))]) %>% rowSums()

write.table(all_abu_frame, alldepth, sep="\t", quote=F, row.names=F)