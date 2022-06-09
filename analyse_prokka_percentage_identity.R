setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)
require(Biostrings)
data("BLOSUM62")

spike_AA <- readAAStringSet("data/proteins/extracted_proteins/all_spikes.faa")
names(spike_AA)[names(spike_AA) == "NC_045512"] <- "SARS-CoV-2"
names(spike_AA)[names(spike_AA) == "NC_004718"] <- "SARS-CoV-1"
names(spike_AA)[names(spike_AA) == "MZ937003"] <- "BANAL-236"
names(spike_AA)[names(spike_AA) == "NC_019843"] <- "MERS"
morsels <- foreach(i = seq(length(spike_AA))) %do% {
  crumbs <- foreach (j = seq(length(spike_AA))) %do% {
    aln <- pairwiseAlignment(spike_AA[[i]], spike_AA[[j]], 
                             substitutionMatrix = BLOSUM62, 
                             gapOpening = 0, gapExtension = -5)
    spike_aa_ident <- pid(aln, type = "PID1")
    accession1 <- names(spike_AA)[i]
    accession2 <- names(spike_AA)[j]
    
    tibble(accession1 = accession1,
           accession2 = accession2,
           pident = spike_aa_ident)
  }
  
  bind_rows(crumbs)
  
}

final <- bind_rows(morsels)

all_viruses <- rev(c("RFGB01", "RFGB02", 
               "RHGB02", "RHGB03",
               "RHGB04", "BANAL-236",
               "SARS-CoV-1", "SARS-CoV-2", 
               "PAGB01", "MERS"))

hm_df <- final %>% 
  filter(accession1 %in% all_viruses,
         accession2 %in% all_viruses)
hm_df %>%
  mutate(accession1 = factor(accession1, all_viruses),
         accession2 = factor(accession2, rev(all_viruses))) %>%
  ggplot(aes(x = accession1, y = accession2, fill = pident)) +
    geom_tile() +
    geom_text(aes(label = round(pident, 1))) +
    scale_fill_gradient(low = "dodgerblue4", high = "paleturquoise2") +
    theme_bw() +
    theme(legend.position = "top",
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(fill = "Amino acid similarity (%; incl. gaps)")

ggsave("results/pident_hm.png", dpi = 600, height = 6, width = 6)

spike_AA[all_viruses]
