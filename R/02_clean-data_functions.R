p_load(tidyverse, phyloseq)

# List of Contaminated Samples i.e. insect falling directly into the cup during eDNA collection
contaminated_samples <- c(
  "F_P02_L1", "F_P02_L2", "E_T03_L1", "E_T03_L2",
  "A_P01_L1", "A_P01_L2", "F_T02_L1", "F_T02_L2"
)

# Function to Trim Sample name
getSampleID <- function(x) {
  x %>%
    str_remove("^Neo_20250321_") %>%
    str_remove("_R1.fastq.gz$")
}

# Sets original sample ID as rownames
trimSampleID <- function(x) {
  x %>%
    mutate(TrimmedSamples = getSampleID(Sample)) %>%
    select(-Sample) %>%
    select(TrimmedSamples, everything()) %>% #<- making sample names first column
    column_to_rownames(var = "TrimmedSamples")
}

# Arthropods with identified class, removing human and fungi samples
arthroFilter <- function(x) {
  x %>%
    as("data.frame") %>% # Coerce into data frame to be accessed by filter f(x)
    filter(phylum == "Arthropoda" & class != "") %>%
    as("matrix") # Coerce into matrix to preserve row and col names in phyloseq
}


clean_phyloseq <- function(x, contam) {
  keep <- setdiff(sample_names(x), contam)
  decontam <- prune_samples(keep, x)
  
  y <- prune_samples(sample_sums(decontam) > 0, decontam)
  y1 <- prune_taxa(taxa_sums(y) > 0, y)
  
  y1 %>%
    subset_samples(Site != "Blank")
}