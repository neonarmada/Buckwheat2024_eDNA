library(sf)
library(ggplot2)
library(tidyverse)
library(lme4)
library(MuMIn)
library(ggsignif)

buckwheat <- read.csv("data/CollectedData_v4.csv") %>%
  drop_na(DevRate)
head(buckwheat)

## Mature Seed Set LM & ANOVA
DevRate_aov <- lm(DevRate ~ ForestArea + Variety + Morph + HerbivRate, data = buckwheat)
summary(DevRate_aov)
DevRate_anova <- anova(DevRate_aov)
DevRate_anova

## Model Selection
DevRate_aov2 <- lm(DevRate ~ ForestArea*Variety*Morph*HerbivRate, data = buckwheat)
m1<-dredge(DevRate_aov2,rank="AICc")

## Mature Seed Set LM & ANOVA (multi-level interactions)
DevRate_aov2.1 <- lm(DevRate ~ ForestArea*Variety*Morph, data = buckwheat)
summary(DevRate_aov2.1)
DevRate_anova1 <- anova(DevRate_aov2.1)
DevRate_anova1

## SeedSet LM & ANOVA
seedset_aov <- lm(SeedSet ~ ForestArea+Variety+Morph+HerbivRate, data = buckwheat)
summary(seedset_aov)
seedset_anova <- anova(seedset_aov)
seedset_anova

## SeedSet LM & ANOVA (multi-level interactions)
seedset_aov1 <- lm(SeedSet ~ ForestArea*Morph*Variety, data = buckwheat)
summary(seedset_aov1)
seedset_anova1 <- anova(seedset_aov1)
seedset_anova1

## Herbivory (multi-level interactions)
herbiv_aov <- lm(HerbivRate ~ ForestArea * Morph * Variety, data = buckwheat)
summary(herbiv_aov)
herbrate_anova <- anova(herbiv_aov)
herbrate_anova

## Herbivory (multi-level interactions)
herbiv_aov1 <- lm(HerbivRate ~ ForestArea * Morph * Variety, data = buckwheat)
summary(herbiv_aov1)
herbrate_anova1 <- anova(herbiv_aov1)
herbrate_anova1

## Site Yield LM & ANOVA
siteyield_aov <- lm(SiteYield ~ ForestArea + Height + Morph + HerbivRate, data = buckwheat)
summary(siteyield_aov)
anova(siteyield_aov)

# Wrtite csv
write.csv(DevRate_anova, "output/DevRate_AOV1.csv")
write.csv(DevRate_anova1, "C:/Users/M9NRC/Desktop/Schoolwork/Statistics/DevRate_multi.csv")
write.csv(seedset_anova, "C:/Users/M9NRC/Desktop/Schoolwork/Statistics/SeedSet_AOV1.csv")
write.csv(seedset_anova1, "C:/Users/M9NRC/Desktop/Schoolwork/Statistics/SeedSet_multi.csv")
write.csv(herbrate_anova, "C:/Users/M9NRC/Desktop/Schoolwork/Statistics/HerbrateAOV1.csv")
write.csv(herbrate_anova1, "C:/Users/M9NRC/Desktop/Schoolwork/Statistics/HerbrateAOV_multi.csv")