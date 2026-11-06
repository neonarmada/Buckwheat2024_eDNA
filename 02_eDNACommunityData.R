#R Script v.11
# Library Call =====
library(sf)
library(ggplot2)
library(tidyverse)
library(lme4)
library(MuMIn)
library(ggsignif)
library(vegan)
library(phyloseq)
library(ggsci)
library(RColorBrewer)
library(fuzzyjoin)

# Import to Phyloseq ===========

# Community Read Data
set.seed(11) # set seed to get reproducible ordination values
community_reads <- read.csv("data/Neo2024_edna_community.csv")

# Filtering Blanks and Samples
blank_reads <- community_reads %>%
  filter(str_detect(Sample, "Negative") | str_detect(Sample, "Blank"))
pull(blank_reads, "Sample")
sample_reads <- community_reads %>%
  filter(!str_detect(Sample, "Negative") & !str_detect(Sample, "Blank"))

# Separating by Site, FlowerMorph, and Leaf/Flower
newcols <- sample_reads %>%
  separate(Sample, sep = "_", into = c("1", "2", "Site", "Morph", "Type")) %>%
  select(Site, Morph, Type) %>%
  # changing letters into words for the categories
  mutate(Morph = case_when(str_starts(Morph, "P") ~ "Pin", str_starts(Morph,"T")
                           ~ "Thrum", TRUE ~ Morph)) %>%
  mutate(Type = case_when(str_starts(Type, "F") ~ "Flower", str_starts(Type,"L")
                          ~ "Leaf", TRUE ~ Type)) # Keep original value if doesn't match

# Adding Category columns back into table
sample_reads2 <- sample_reads %>%
  add_column(newcols, .after = "Sample")

# Blank Reads Category
head(blank_reads)
blank_reads2 <- blank_reads %>%
  mutate(Type = case_when(str_detect(Sample, "Negative") ~ "Control", 
                          str_detect(Sample, "Blank") ~ "Blank", TRUE ~ NA), .after = Sample)

# Recombine Blank + Samples
community_reads2 <- bind_rows(sample_reads2, blank_reads2)

# Trimming Sample Names 
community_reads2$TrimmedSamples <- str_remove(community_reads2$Sample, "^Neo_20250321_") %>%
  str_remove("_R1.fastq.gz$")
community_reads2 <- community_reads2 %>%
  select(-Sample)

# Group and Arranging samples by Site
community_reads3 <- community_reads2 %>%
  select(TrimmedSamples, everything()) %>% #<- making sample names first column
  group_by(Site) %>%
  column_to_rownames(var = "TrimmedSamples")

# Import to OTU Table
community_otu <- community_reads3 %>% 
  select(-Site, -Morph, -Type)
otu <- otu_table(community_otu, taxa_are_rows = FALSE)

# Import to Taxonomy Table
edna_taxon <- read.csv("data/Neo2024_edna_taxon_v2.csv") %>%
  column_to_rownames(var = "query") %>%
  select(-No) %>%
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
buckwheat_filtered <- prune_samples(sample_sums(buckwheat_raw)>0, buckwheat_raw) # removing samples with low total reads
buckwheat_filtered1 <- subset_samples(buckwheat_filtered, Site != "Blank") #removing blanks
buckwheat_filtered2 <- subset_samples(buckwheat_filtered1, !sample_names(buckwheat_filtered1) %in% #Filtering Out Suspected Cases of Contamination i.e. insect falling directly into the cup
                                        c("F_P02_L1", "F_P02_L2","E_T03_L1", "E_T03_L2",
                                          "A_P01_L1", "A_P01_L2", "F_T02_L1", "F_T02_L2"))
# Ordination========================
# NMDS Ordination
nmds_filt <- ordinate(buckwheat_filtered2, "NMDS", "bray", autotransform = F)
plot(nmds_filt)
nmds_result_filt <- as.data.frame(scores(nmds_filt,display=c("sites")))
spscrs_filt <- as.data.frame(scores(nmds_filt,display=c("species")))

# Join NMDS to Sample Data
nmds_result_filt2 <- nmds_result_filt %>%
  rownames_to_column(var = "Samples")
dat2 <- sample_data %>%
  rownames_to_column(var = "Samples")
nmds_result_comb <- left_join(nmds_result_filt2, dat2, by = "Samples")

#Convex Hull (adding the polygon shape to NMDS)
nmds_hull <- nmds_result_comb %>%
  group_by(Site) %>%
  slice(chull(NMDS1,NMDS2))

# NMDS Plot
nmds_plot_filt <- ggplot(nmds_result_comb, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = Site), shape =21, size = 3, stroke = 0.7)+
  geom_polygon(data = nmds_hull, alpha = 0.2, aes(fill = Site,colour = Site)) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2")+
  theme_bw() +
  theme_bw(base_size = 28) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
nmds_plot_filt
ggsave("output/edna_filt_nmds.png",
       plot = nmds_plot_filt, width = 2400 , height = 1600, units = "px", dpi = 300)

# Species Scores
spscrs_filt
spscrs_plot_filt <- ggplot(spscrs_filt, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = .5, stroke = 1) +
  theme_bw() +
  theme_bw(base_size = 28) +
  theme(legend.position = "none", axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
spscrs_plot_filt
ggsave("output/edna_filt_spscores.png",
       plot = spscrs_plot_filt, width = 2400, height = 1600, units = "px", dpi = 300)

# Fix Anova!
# Adonis ==================================

adonis_otu <- otu_table(buckwheat_filtered2)
adonis_dat <- sample_data(buckwheat_filtered2)

# Permanova Analysis
permanova2 <- adonis2(adonis_otu ~ Site, adonis_dat, perm=9999, by="terms")
permanova2

write.csv(permanova2, "output/eDNA_adonis3.csv")

# Yield/Herbiv x Community ================
### Seed Set Linear Model

veg_nmds_result2 <- nmds_result_filt %>%
  rownames_to_column(var = "Samples")
nmds_veg_comb <- left_join(veg_nmds_result2, dat2, by = "Samples")
lm_df <- nmds_veg_comb #nmds_result_comb

seedset_lm <- lm(SeedSet ~ NMDS1 + NMDS2, data = lm_df)
summary(seedset_lm)
seedset_anova <- anova(seedset_lm)
seedset_anova # Nothing

### Dev Rate Linear Model
devrate_lm <- lm(DevRate ~ NMDS1 + NMDS2, data = lm_df)
summary(devrate_lm) # NMDS2 Negatively Assoc. with Dev Rate
devrate_anova <- anova(devrate_lm)
devrate_anova # NMDS2 Negatively Assoc. with Dev Rate !!!!!

devrate_NMDS2 <- ggplot(lm_df, aes(x = NMDS2, y = DevRate)) +
  geom_point(size = .5, stroke = 1)

### Herbiv Linear Model
herbiv_lm <- lm(HerbivRate ~ NMDS1 + NMDS2, data = lm_df)
summary(herbiv_lm)
herbiv_anova <- anova(herbiv_lm)
herbiv_anova # Nothing

# Misc ======
#final phyloseq ordin_filt
#final otu version veg_otu3
tax_filt <- tax_table(ordin_filt)
