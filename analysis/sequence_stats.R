#Sequence_stats

#Libraries
library(tidyverse)
library(data.table)
library(kableExtra)

#Files
sample_reads <- read_tsv("../sequencing_reports/J.Pinhassi_21_02_sample_info.txt") %>%
  rename(sample = 1)

sample_ID <- read_tsv("../data/J.Pinhassi_21_02_sample_info.txt") %>%
  dplyr::select(1:2) %>%
  dplyr::rename(sample = 1, sample_name = 2) %>%
  mutate(sample_name = str_replace(sample_name, "[0-9]*:","")) %>%
  separate(sample_name, c("treatment","timepoint"), sep = ",") %>%
  mutate(treatment = gsub("NT","TN", treatment)) %>%
  mutate(sample_name = paste(timepoint, treatment, sep = "_")) %>%
  separate(treatment, c("treatment","replicate"), sep = -1) %>%
  mutate(tre_rep = gsub("_","", sample_name)) %>%
  mutate(tre_rep = toupper(tre_rep)) %>%
  mutate(day = case_when( timepoint == "t3" ~ 10,
                          TRUE ~ 17
  )
  )

overall_stats <- read_tsv("../data/overall_stats.tsv") %>%
  dplyr::select(-2,-3)  %>%
  mutate(sample = str_extract(sample, "P[0-9]*_[0-9]*"))

#Reading in taxonomic annotation of reads
taxonomy <- read_tsv("../data/eukulele_phylodb.tsv") %>%
  select(-species) %>%
  rename(species = "genus", genus = "family", family = "order", order = "class", class = "phylum", phylum = "kingdom")

#Count file for all orfs, tpm is calculated per sample
bbmap <- fread("../data/bbmap_counts.tsv.gz", sep = "\t") %>%
  filter(count > 0)  %>%                               # removing 0 counts to reduce size of table
  mutate(Geneid = str_replace(Geneid, "[0-9]*_",""),   # Joining in chr and geneid to match format of eggnog naming.
         orf = paste(Chr,Geneid, sep = "_"),
         sample = str_extract(sample, "P[0-9]*_[0-9]*")) %>% # Removing extra numbers in sample to match sample_ID format
  dplyr::select(-Geneid,-Chr) %>%
  dplyr::select(orf,Start,End,Strand,Length,sample,count,tpm) 

#Putting it all together and writing the file
sample_ID %>%
  inner_join(sample_reads, by = "sample") %>%
  select(sample, treatment, timepoint, replicate, Mreads) %>%
  inner_join(overall_stats, by = "sample") %>%
  inner_join(
    bbmap %>%
      inner_join(taxonomy, by = "orf") %>%
      group_by(sample, domain) %>%
      summarise(counts = sum(count)) %>%
      ungroup() %>%
      filter(domain !="NA") %>%
      spread(domain, counts, fill = 0),
    by = "sample"
  ) %>%
  write_tsv(.,"../results/sequencing_stats.tsv")
