p_load(tidyverse, phyloseq)
set.seed(123)

# Site Letter to Code
site_key <- tribble(
  ~from, ~to,
  "A", "SJ",
  "B", "SW",
  "C", "UN",
  "D", "SG",
  "E", "UB",
  "F", "TH",
) %>%
  mutate(to = as.character(to))

# Function: Replacing vague B4 site names
replace_site <- function(x) {
  x %>%
    mutate(Site = replace_values(Site,
                                        from = site_key$from,
                                        to = site_key$to
    ))
}

# 01_YieldData ====
# Calculating flower estimates based on LM
combine_lm <- function(data, lm) {
  # Extract p values from anova
  p_vals <- lm %>%
    anova() %>%
    as.data.frame() %>%
    mutate(PrF = `Pr(>F)`) %>%
    select(PrF) %>%
    arownames_to_column(var = "coefficient")

  # Extract coefficients from lm
  coeff <- as.data.frame(summary(lm)$coefficients) %>%
    rownames_to_column(var = "coefficient") %>%
    full_join(p_values, by = c("coefficient")) %>%
    filter_out(is.na(Estimate)) %>%
    filter_out(PrF >= 0.05) # filter only values less than or eq to 0.05, this doesnt work if intercept is not signif

  # Multiply measured traits by LM estimates
  data_lm <- data %>%
    mutate(
      across(any_of(coeff$coefficient), # check all the coefficients available in coeff
             function(x) {
               x * coeff$Estimate[match(cur_column(), coeff$coefficient)] # multiply values by coefficients estimates
               }, 
             .names = "{.col}_est") # regular naming expression "(colname)_est"
      )
  
  # Sum with intercept to get flower count
  data_est <- data_lm %>%
    mutate(estimate = coeff$Estimate[match("(Intercept)", coeff$coefficient)] +
             + rowSums(pick(ends_with("_est")), na.rm = TRUE)
      )
  
  # Coalesce no. of flowers into one column (prioritising original flower count)
  data_est %>%
    mutate(mixed = round(coalesce(data_est$Flower, data_est$estimate)))
}

plant_yield <- function(plant, seed, flowr) {
  herb <- plant_traits %>%
    mutate(HerbivRate = DmgLeaf/TotalLeaf)
  
  seed_data %>%
    select(PlantID, everything())%>%
    full_join(inflo_sum, "PlantID") %>%
    mutate(DevRate = Dev/flower) %>%
    mutate(SeedSet = (Dev+Undev)/flower) %>%
    mutate(AvgMass = TotalDevMass/Dev) %>%
    full_join(herb, "PlantID")
  
}


# 02_eDNACommunityData ====
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
