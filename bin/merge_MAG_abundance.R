library(data.table)
library(tidyverse)

mag_dt <- rbindlist(lapply(list.files(pattern= "_abundance_table.tbl"), fread))

mag_dt$SampleID <- gsub("_cleanbin_[0-9]+","",mag_dt$bin)


mag_abudance_wide <- pivot_wider(mag_dt[,-"bin"],names_from = "SampleID", values_from = "scaled_weighted_contigs_TPM_per_bin", values_fill = 0)

fwrite(mag_abudance_wide, file="merged_MAG_tpm_abundance.tbl")

taxsplit = strsplit(mag_abudance_wide$taxa, split=";", fixed=TRUE)
taxunitmatrix = matrix(NA, ncol=max(sapply(taxsplit, length)), nrow=length(taxsplit))
colnames(taxunitmatrix) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxunitmatrix)]
rownames(taxunitmatrix) = mag_abudance_wide$taxa
for (i in 1:nrow(taxunitmatrix)){
    taxunitmatrix[i, 1:length(taxsplit[[i]])] <- taxsplit[[i]]
}
taxunitmatrix = gsub("[a-z]__", "", taxunitmatrix)

phylum_dt <- data.frame(Phylum=taxunitmatrix[,"Phylum"], mag_abudance_wide) %>% group_by(Phylum) %>%
    summarise(across(-c(taxa), sum)) %>%  pivot_longer(cols = -c("Phylum"), names_to = "SampleID",values_to = "Abundance")


phylum_abudance_plot <- ggplot(phylum_dt,aes(fill=Phylum,x=SampleID,y=Abundance))+
    geom_bar(position="stack",stat="identity")+
    theme(axis.text.x = element_text(angle=90))

ggsave(phylum_abudance_plot,file="phylum_rel_abudance_plot.png")