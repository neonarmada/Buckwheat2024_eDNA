# R Script v.12
# Library Call =====
p_load(ggplot2, tidyverse, vegan, phyloseq, lme4, MuMIn, styler, RColorBrewer, viridis)
library(sf)
library(ggsignif)
library(ggsci)

# Import to Phyloseq ===========
set.seed(123)
source("R/01_load-data_functions.R")
source("R/02_clean-data_functions.R")
source("R/04_generate-figs_functions.R")

# Importing OTU Table
community_otu <- trimSampleID(community_reads)
otu <- otu_table(community_otu, taxa_are_rows = FALSE)

# Importing TAX Table
community_tax <- arthroFilter(edna_taxon)
tax <- tax_table(community_tax)

# Importing DAT Table
dat <- sample_data(community_info)

# Combine into phyloseq object
buckwheat_raw <- phyloseq(otu, tax, dat)
buckwheat_raw

# Phyloseq Pre-Processing
buckwheat_filtered <- clean_phyloseq(buckwheat_raw, contaminated_samples) 

# Exploratory Analysis =====
## NMDS Ordination ======
nmds_filt <- ordinate(buckwheat_filtered, "NMDS", "bray", autotransform = F)
nmds_scrs <- scores(nmds_filt, display = c("sites")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Samples")
spscrs_filt <- scores(nmds_filt, display = c("species")) %>%
  as.data.frame()

dat2 <- community_info %>%
  rownames_to_column(var = "Samples")
nmds_scrs_comb <- inner_join(nmds_scrs, dat2, by = "Samples")

## ANOVA/ANOSIM ========
adonis_otu <- otu_table(buckwheat_filtered) %>%
  data.frame()
adonis_dat <- sample_data(buckwheat_filtered) %>%
  data.frame()

# ANOVA from community data matrix
otu_anova <- adonis2(adonis_otu ~ ForestArea + Type + Morph + Variety, data = adonis_dat, perm = 9999, by = "terms")
otu_anova #variety not so accurate maybe?

# ANOSIM from community data matrix
anosim(adonis_otu, adonis_dat$ForestArea, distance = "bray", perm = 9999)
anosim(adonis_otu, adonis_dat$Type, distance = "bray", perm = 9999)
anosim(adonis_otu, adonis_dat$Morph, distance = "bray", perm = 9999)
anosim(adonis_otu, adonis_dat$Variety, distance = "bray", perm = 9999)

# Linear Model for NMDS1&2 (I dont know if this is correct or not)
lm(NMDS1 ~ ForestArea + Type + Morph + Variety, data = nmds_scrs_comb) %>%
  anova()
lm(NMDS2 ~ ForestArea + Type + Morph + Variety, data = nmds_scrs_comb) %>%
  anova() #NMDS2 significantly influences ForestArea & Type

## Yield/Herbiv x Community ================
### Seed Set Linear Model
seedset_lm <- lm(SeedSet ~ NMDS1 + NMDS2, data = nmds_scrs_comb)
summary(seedset_lm)
seedset_anova <- anova(seedset_lm)
seedset_anova # NMDS1 assoc w. seedset sometimes depending on set.seed() value

### Dev Rate Linear Model
devrate_lm <- lm(DevRate ~ NMDS1 + NMDS2, data = nmds_scrs_comb)
summary(devrate_lm) 
devrate_anova <- anova(devrate_lm)
devrate_anova # NMDS2 Negatively Assoc. with Dev Rate !!!!!

### Herbiv Linear Model
herbiv_lm <- lm(HerbivRate ~ NMDS1 + NMDS2, data = nmds_scrs_comb)
summary(herbiv_lm)
herbiv_anova <- anova(herbiv_lm)
herbiv_anova # NMDS2 Positively Assoc with herbiv Rate !!!!

# Data Visualisation ====
# Convex Hull (adding the polygon shape to NMDS)
nmds_hull <- nmds_scrs_comb %>%
  group_by(Site) %>%
  slice(chull(NMDS1, NMDS2))

# NMDS Plot
nmds_plot <- ggplot(nmds_scrs_comb, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(fill = Site), shape = 21, size = 3, stroke = 0.7) +
  geom_polygon(data = nmds_hull, alpha = 0.2, aes(fill = Site, colour = Site)) +
  scale_fill_viridis(end=0.9,direction=-1,option="magma",discrete = T) +
  scale_colour_viridis(end=0.9,direction=-1,option="magma",discrete = T) +
  theme_bw(base_size = 28) +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black")
  )
nmds_plot

devrate_NMDS2 <- ggplot(nmds_scrs_comb, aes(x = NMDS2, y = DevRate)) +
  geom_point(size = .5, stroke = 1)
devrate_NMDS2

herbiv_NMDS2 <- ggplot(nmds_scrs_comb, aes(x = NMDS2, y = HerbivRate)) +
  geom_point(size = .5, stroke = 1)
herbiv_NMDS2