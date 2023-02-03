library(phyloseq)
library(tidyverse)
library(speedyseq)
allBootstrappedData <- function(tl, rarefy.depth=1080){
  # Check rarefaction depth

  rarefied <- rarefy_even_depth(tl, sample.size=rarefy.depth, verbose=FALSE)
  psm<-speedyseq::psmelt(rarefied)

}

lan_diff<-readRDS('/Users/kalen/Documents/YNP_Project/16S/SC_Overall/Overall_BA_GCN.RDS') #Update to correct path
lan_diff<-prune_samples(sample_sums(lan_diff)>1079, lan_diff) #make sure this value matches your rarefaction depth or that it is at least less than it.
lan_diff<-phyloseq::filter_taxa(lan_diff, function(x) max(x)>1, TRUE)
dc<-do.call("bind_rows", replicate(65, allBootstrappedData(lan_diff), simplify = FALSE))%>% group_by(OTU, Sample) %>% dplyr::mutate( Abundance = mean(Abundance), stDev= sd(Abundance)) %>% distinct()%>%pivot_wider(id_cols='Sample', names_from='OTU', values_from='Abundance')
#Update the above like to contain the number of iterations to perform under the replicate function.

#dc2<-dc%>% group_by(OTU, Sample)%>% mutate(Abundance=mean(Abundance))


#dc3<-dc2 %>% distinct()
#dc4<-dc3%>%pivot_wider(id_cols='Sample', names_from='OTU', values_from='Abundance')
dc4<-column_to_rownames(dc, var = "Sample")
boot<-phyloseq(otu_table(dc4, taxa_are_rows=FALSE), sample_data(sample_data(lan_diff)), tax_table(tax_table(lan_diff)))
#boot<-phylosmith::relative_abundance(boot)
#boot<- transform_sample_counts(boot, function(x) x*100 )
#boot<-filter_taxa(boot, function(x) sum(x > 0.3) > (0.05*length(x)), TRUE)
saveRDS(boot,'/Users/kalen/Documents/YNP_Project/16S/SC_Overall/SC_multirare.RDS') #Update to correct path
