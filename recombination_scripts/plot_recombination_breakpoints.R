setwd("../Desktop/git_repos/bat-CoVs/")
require(tidyverse)
require(data.table)
require(foreach)

# Gene annotations
gene_df <- fread("results/recombination_out/rdp_out_masked/NC_025217_genes.gff") %>% 
  select(V2, type = V3, start = V4, end = V5, info = V9) %>%
  filter(type %in% c("gene", "three_prime_UTR", "five_prime_UTR")) %>%
  separate(info, into = c(rep(NA, 2), "gene_name"), sep = ";")

parsed_gene <- gene_df %>%
  mutate(gene_name = case_when(type == "three_prime_UTR" ~ "3'UTR",
                               type == "five_prime_UTR" ~ "5'UTR",
                               TRUE ~ gene_name)) %>%
  mutate(gene_name = gsub("Name=", "", gene_name)) %>%
  # Merge Spike CDS
  mutate(start = ifelse(gene_name == "S", 22058, start),
         end = ifelse(gene_name == "S", 27516, end)) %>%
  filter(gene_name != "NA39_gp2")

# Recombination events
df <- fread("results/recombination_out/rdp_out_masked/Sarbecovirus.masked.parsed.csv", header = T)
colnames(df) <- c("event_no", "event_no_in_rdp", "aln_start", 
                  "aln_end", "seq_start", "seq_end", "rel_start",
                  "rel_end", "sequences", "minor_parents", "major_parents",
                  "rdp", "geneconv", "bootscan", "maxchi",
                  "chimaera", "sisscan", "phylpro", "lard",
                  "threeseq")

parsed <- df %>% 
  select(event_no, sequences, rel_start, rel_end, rdp, geneconv, 
         bootscan, maxchi, chimaera, sisscan, phylpro, lard, threeseq)

# Get only events supported by all methods
bp <- parsed %>% 
  filter(rel_start != "") %>%
  pivot_longer(c("rdp", "geneconv", "bootscan", "maxchi",
                 "chimaera", "sisscan", "phylpro", "lard",
                 "threeseq"), names_to = "method", values_to = "pval") %>%
  group_by(event_no, rel_start, rel_end) %>%
  summarise(n_sig = sum(pval != "NS")) %>%
  filter(n_sig >= 6)

bp_morsels <- foreach(event = unique(bp$event_no), .combine = "c") %do% {
  sstart <- bp[bp$event_no == event, ]$rel_start
  sstart <- as.numeric(gsub("[[:punct:]]", "", sstart))
  send <- bp[bp$event_no == event, ]$rel_end
  send <- as.numeric(gsub("[[:punct:]]", "", send))
  seq(sstart, send)
}

# Get gene annotation plot
gene_annot <- parsed_gene %>%
  mutate(gene_name = factor(gene_name, parsed_gene$gene_name)) %>%
  ggplot() +
  geom_rect(aes(xmin = start, xmax = end, 
                ymin = 0, ymax = 0.5,
                fill = gene_name),
            color = "black",
            linewidth = 0.5) +
  scale_x_continuous(breaks = c(1, 5000, 10000, 15000, 20000, 25000, 30000)) +
  geom_text(aes(x = start + (end - start) / 2,
                y = 0.25, 
                label = gene_name), angle = 90) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

# Get dashed lines for gene breaks
long_pos <- parsed_gene %>%
  select(start, end) %>%
  pivot_longer(everything(), names_to = "type", values_to = "pos")

dist2 <- tibble(pos = bp_morsels) %>%
  group_by(pos) %>%
  summarise(freq = n()) %>%
  ggplot(aes(x = pos, y = freq, fill = freq)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = c(1, 5000, 10000, 15000, 20000, 25000, 30000),
                     position = "top") +
  geom_vline(xintercept = long_pos$pos, 
             lty = "dashed",
             size = 0.2) +

  scale_fill_gradient2(low = "blue",
                       mid = "khaki1",
                       high = "red",
                       midpoint = 20) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  labs(x = "Pos. on NC_025217 (root)", y = "Num. of events")

ggpubr::ggarrange(dist2, gene_annot, 
                  nrow = 2,
                  heights = c(7, 1),
                  align = "v")

ggsave("results/recombination_out/rdp4_breakpoint_distribution.pdf", dpi = 600, 
       width = 10.5, height = 6)

# Get novel genome recombinants
filt_events <- parsed %>% 
  filter(rel_start != "") %>%
  pivot_longer(c("rdp", "geneconv", "bootscan", "maxchi",
                 "chimaera", "sisscan", "phylpro", "lard",
                 "threeseq"), names_to = "method", values_to = "pval") %>%
  group_by(event_no, rel_start, rel_end) %>%
  summarise(n_sig = sum(pval != "NS"))

novel_events <- df %>% 
  filter(grepl("RFGB01|RFGB02|RHGB01|RHGB02|MW719567", major_parents) |
           grepl("RFGB01|RFGB02|RHGB01|RHGB02|MW719567", minor_parents) |
           grepl("RFGB01|RFGB02|RHGB01|RHGB02|MW719567", sequences))


novel_df <- filt_events %>% 
  filter(event_no %in% novel_events$event_no) %>%
  mutate(rel_end = as.numeric(gsub("\\*", "", rel_end)),
         rel_start = as.numeric(gsub("\\*", "", rel_start))) %>%
  arrange(rel_end) %>%
  add_column(offs = c(0.5, 0, 0.5, 0, 0.5, 0))
  
novel_plt <- novel_df %>%
  ggplot() +
  geom_rect(aes(xmin = rel_start, xmax = rel_end, 
                ymin = offs + 0, ymax = offs + 0.5,
                fill = n_sig),
            color = "black",
            linewidth = 0.5) +
  scale_x_continuous(breaks = c(1, 5000, 10000, 15000, 20000, 25000, 30000, 35000)) +
  theme(legend.position = "bottom", 
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) 
novel_df

ggpubr::ggarrange(dist2, gene_annot, novel_plt, 
                  nrow = 3,
                  heights = c(7, 1, 2),
                  align = "v")

ggsave("results/recombination_out/novel_breakpoints.pdf", dpi = 600, 
       width = 10.5, height = 8)

