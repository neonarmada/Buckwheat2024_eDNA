# R Script v.12
# Library Call =====
p_load(ggplot2, tidyverse, vegan, phyloseq, lme4, MuMIn)
library(sf)
library(ggsignif)
library(ggsci)
library(RColorBrewer)

# Import to Phyloseq ===========

# Community Read Data
set.seed(123) # set seed to get reproducible ordination values
community_reads <- read_delim("data/raw/Neo2024_edna_community.tsv", delim = " ") %>%
  rename(Sample = query)

# Function to Trim Sample name
getSampleID <- function(x) {
  x %>%
    str_remove("^Neo_20250321_") %>%
    str_remove("_R1.fastq.gz$")
}

# Preparing OTU Table
community_otu <- community_reads %>%
  mutate(TrimmedSamples = getSampleID(Sample)) %>%
  select(-Sample) %>%
  select(TrimmedSamples, everything()) %>% #<- making sample names first column
  column_to_rownames(var = "TrimmedSamples")

otu <- otu_table(community_otu, taxa_are_rows = FALSE)

# Import to Taxonomy Table
edna_taxon <- read_tsv("data/raw/taxonomy_qc_Neo2024_arthropod.tsv") %>%
  column_to_rownames(var = "query") %>%
  filter(phylum == "Arthropoda" & class != "") # Remove human and fungi reads, only arthropods with identified class

tax <- tax_table(edna_taxon)
rownames(tax) <- rownames(edna_taxon)
colnames(tax) <- colnames(edna_taxon)

# Import Sample Information
sample_data <- read.csv("data/CollectedData_v5.csv") %>%
  column_to_rownames(var = "ID")
dat <- sample_data(sample_data) # <- Imports as Phyloseq sample data table

# Combine into phyloseq object
buckwheat_raw <- phyloseq(otu, tax, dat)
buckwheat_raw

# Pre-Processing
buckwheat_filtered <- prune_samples(sample_sums(buckwheat_raw) > 0, buckwheat_raw) # removing samples with low total reads
buckwheat_filtered1 <- subset_samples(buckwheat_filtered, Site != "Blank") # removing blanks
buckwheat_filtered2 <- subset_samples(buckwheat_filtered1, !sample_names(buckwheat_filtered1) %in% # Filtering Out Suspected Cases of Contamination i.e. insect falling directly into the cup
  c(
    "F_P02_L1", "F_P02_L2", "E_T03_L1", "E_T03_L2",
    "A_P01_L1", "A_P01_L2", "F_T02_L1", "F_T02_L2"
  ))
# Ordination========================
# NMDS Ordination
nmds_filt <- ordinate(buckwheat_filtered2, "NMDS", "bray", autotransform = F)
plot(nmds_filt)
nmds_result_filt <- as.data.frame(scores(nmds_filt, display = c("sites")))
spscrs_filt <- as.data.frame(scores(nmds_filt, display = c("species")))

# Join NMDS to Sample Data
nmds_result_filt2 <- nmds_result_filt %>%
  rownames_to_column(var = "Samples")
dat2 <- sample_data %>%
  rownames_to_column(var = "Samples")
nmds_result_comb <- left_join(nmds_result_filt2, dat2, by = "Samples")

# Convex Hull (adding the polygon shape to NMDS)
nmds_hull <- nmds_result_comb %>%
  group_by(Site) %>%
  slice(chull(NMDS1, NMDS2))

# NMDS Plot
nmds_plot_filt <- ggplot(nmds_result_comb, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = Site), shape = 21, size = 3, stroke = 0.7) +
  geom_polygon(data = nmds_hull, alpha = 0.2, aes(fill = Site, colour = Site)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() +
  theme_bw(base_size = 28) +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  )
nmds_plot_filt
ggsave("output/edna_filt_nmds.png",
  plot = nmds_plot_filt, width = 2400, height = 1600, units = "px", dpi = 300
)

# Species Scores
spscrs_filt
spscrs_plot_filt <- ggplot(spscrs_filt, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = .5, stroke = 1) +
  theme_bw() +
  theme_bw(base_size = 28) +
  theme(
    legend.position = "none", axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  )
spscrs_plot_filt
ggsave("output/edna_filt_spscores.png",
  plot = spscrs_plot_filt, width = 2400, height = 1600, units = "px", dpi = 300
)

# Adonis ==================================
adonis_otu <- otu_table(buckwheat_filtered2) %>%
  data.frame()
adonis_dat <- sample_data(buckwheat_filtered2) %>%
  data.frame()

# Permanova Analysis
permanova2 <- adonis2(adonis_otu ~ ForestArea + Type + Morph + Variety, data = adonis_dat, perm = 9999, by = "terms")
permanova2 

write.csv(permanova2, "output/eDNA_adonis3.csv")

# Yield/Herbiv x Community ================
### Seed Set Linear Model

veg_nmds_result2 <- nmds_result_filt %>%
  rownames_to_column(var = "Samples")
nmds_veg_comb <- left_join(veg_nmds_result2, dat2, by = "Samples")
lm_df <- nmds_veg_comb # nmds_result_comb

seedset_lm <- lm(SeedSet ~ NMDS1 + NMDS2, data = nmds_veg_comb)
summary(seedset_lm)
seedset_anova <- anova(seedset_lm)
seedset_anova # Nothing

### Dev Rate Linear Model
devrate_lm <- lm(DevRate ~ NMDS1 + NMDS2, data = nmds_veg_comb)
summary(devrate_lm) # NMDS2 Negatively Assoc. with Dev Rate
devrate_anova <- anova(devrate_lm)
devrate_anova # NMDS2 Negatively Assoc. with Dev Rate !!!!!

devrate_NMDS2 <- ggplot(nmds_veg_comb, aes(x = NMDS2, y = DevRate)) +
  geom_point(size = .5, stroke = 1)
devrate_NMDS2

### Herbiv Linear Model
herbiv_lm <- lm(HerbivRate ~ NMDS1 + NMDS2, data = nmds_veg_comb)
summary(herbiv_lm)
herbiv_anova <- anova(herbiv_lm)
herbiv_anova # NMDS2 Positively Assoc with herbiv Rate !!!!

herbiv_NMDS2 <- ggplot(nmds_veg_comb, aes(x = NMDS2, y = HerbivRate)) +
  geom_point(size = .5, stroke = 1)
herbiv_NMDS2

# Misc ======
# final phyloseq ordin_filt
# final otu version veg_otu3
tax_filt <- tax_table(ordin_filt)
