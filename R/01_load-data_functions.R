p_load(tidyverse)

# Import ASV Table
community_reads <- read_delim("data/raw/Neo2024_edna_community.tsv", delim = " ") %>%
  rename(Sample = query)

# Import Taxa Table
edna_taxon <- read_tsv("data/raw/taxonomy_qc_Neo2024_arthropod.tsv") %>%
  column_to_rownames(var = "query")

# Import Sample Info
community_info <- read.csv("data/CollectedData_v5.csv") %>%
  column_to_rownames(var = "ID")
