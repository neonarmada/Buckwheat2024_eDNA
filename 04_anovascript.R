library(sf)
library(ggplot2)
library(tidyverse)
library(lme4)
library(MuMIn)
library(ggsignif)
library(broom)
set.seed(123)

## Import Data and Cleanup ====
buckwheat_data <- read.csv("data/CollectedData_v6.csv") %>%
  drop_na(DevRate)
head(buckwheat_data)

### Remove eDNA replicates ====
buckwheat_data$PlantID <- str_sub(buckwheat_data$ID, end=-4) #Remove eDNA suffix
buckwheat <- buckwheat_data %>% #Remove duplicates
  select(PlantID, everything(), -ID, -Type) %>% #remove eDNA ID and type columns and Move PlantID to front
  distinct() #remove duplicated rows


## Mature Seed Set LM & ANOVA (multi-level interactions)
DevRate_LM <- lm(DevRate ~ ForestArea+Variety+Morph, data = buckwheat)
summary(DevRate_LM)
DevRate_anova <- anova(DevRate_LM)
DevRate_anova

## SeedSet LM & ANOVA (multi-level interactions)
seedset_LM <- lm(SeedSet ~ ForestArea+Morph+Variety, data = buckwheat)
summary(seedset_LM)
seedset_anova <- anova(seedset_LM)
seedset_anova

## Herbivory (multi-level interactions)
herbiv_LM <- lm(HerbivRate ~ ForestArea * Morph * Variety, data = buckwheat)
summary(herbiv_LM)
herbrate_anova <- anova(herbiv_LM)
herbrate_anova

## Site Yield LM & ANOVA
siteyield_LM <- lm(SiteYield ~ ForestArea + Variety + HerbivRate + Morph, data = buckwheat)
summary(siteyield_LM)
siteyield_anova <- anova(siteyield_LM)
siteyield_anova

# Write csv
write.csv(DevRate_anova, "output/AOVTables/DevRate_AOV1.csv")
write.csv(seedset_anova, "output/AOVTables/SeedSet_AOV1.csv")
write.csv(herbrate_anova, "output/AOVTables/HerbRate_AOV1.csv")
write.csv(siteyield_anova, "output/AOVTables/SiteYield_AOV1.csv")


