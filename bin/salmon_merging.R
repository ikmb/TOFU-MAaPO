# Based on scripts written by Malte Rühlemann. Adapted by Eike Matthias Wacker 2023

# Script takes all {sample].quant.sf files, calculates the TPM values based on contig lengths and merges them together. If a reference file exists, TPM values will be given for every taxonomic level.

library(dplyr)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
#args <- c("salmon/GloHuGG.GTDBr214.cluster_final_tax.tsv")

if(length(args) >= 1){
  reference = args[1]
  perform_taxing <- TRUE
}else{
  perform_taxing <- FALSE
  reference <- NULL
}

quantfiles <- list.files(path=".", pattern = ".sf")

if(length(quantfiles) == 0) { stop("Did not find any salmon quant output files in this directory")}

for(i in 1:length(quantfiles)){
  sample_name<- quantfiles[i] %>% tools::file_path_sans_ext() %>% tools::file_path_sans_ext()
  ##alternative, if tools library is not available:
  # sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(quantfiles[i]))
  
  if(i == 1){
    left <- fread(paste0("salmon/",quantfiles[i]))
    
    left <- left %>% mutate(Coverage=NumReads*300/EffectiveLength)  %>% 
      mutate(EffectiveLengthCov=ifelse(Coverage>0.1, 1 , 0)*EffectiveLength) %>% 
      select(-NumReads, -Coverage, -EffectiveLengthCov)

    left <- left %>% rename(!!paste0("TPM_",sample_name) := TPM) 
  }
  if(i >= 2){
    otherfile <- fread(paste0("salmon/",quantfiles[i]))
    
    otherfile <- otherfile %>% mutate(Coverage=NumReads*300/EffectiveLength)  %>% 
      mutate(EffectiveLengthCov=ifelse(Coverage>0.1, 1,0)*EffectiveLength) %>% 
      select(-NumReads, -Coverage, -EffectiveLengthCov, -Length, -EffectiveLength) %>% 
      rename(!!paste0("TPM_",sample_name) := TPM) 
    
    left <- left_join(left, otherfile, by="Name")
  }
}

if(!perform_taxing){
  #if no reference file is available, just write the merged tpm file
  fwrite(left, "salmon_merged_TPM.tbl")
  
}else{
  loadedreference <- fread(reference)
  
  colnames(loadedreference) = c("bin","taxa")
  
  taxas = sapply(paste0(loadedreference$bin,";",loadedreference$taxa), function(x){
    x_split=strsplit(x, split=";")[[1]]
    y = x_split[2:8]
    y_low=y[max(grep("__$",y, invert=T))]
    if(max(grep("__$",y, invert=T)) == 6){
      y[7] = paste0(gsub("^g__","s__", y_low),"_sp_", gsub("_cluster95_","_",x_split[1]))
    }
    if(max(grep("__$",y, invert=T)) == 5){ # nolint
      y[6] = paste0("g__",gsub("_cluster95_","_",x_split[1]))
      y[7] = paste0("s__",gsub("_cluster95_","",x_split[1]),"_sp_", gsub("_cluster95_","_",x_split[1]))
    }
    return(gsub(" ","_", c(x_split[1],y)))
  }) %>%
    t %>% data.frame(row.names=NULL, stringsAsFactors=F)
  
  colnames(taxas) = c("bin","kingdom","phylum","class","order","family","genus","species")
  
  tpm_tax <- merge(taxas, left %>% mutate(bin = gsub("_contig_[0-9]+", "", Name)), by = "bin", all.y = T)
  
  complete_tax_table = lapply(c("kingdom","phylum","class","order","family","genus","species","bin"), function(x) return(tpm_tax %>% 
                                                                                                                           select(tax=x, starts_with("TPM")) %>% 
                                                                                                                           group_by(tax) %>% summarise(across(starts_with("TPM"), ~ sum(.x, na.rm = T))) %>% 
                                                                                                                           mutate(level=x) %>% select(tax, level, everything()))) %>% 
    do.call("bind_rows",.)
  
  fwrite(complete_tax_table, "salmon_merged_TPM.tbl")
}
