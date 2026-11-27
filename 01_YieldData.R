p_load(sf, ggplot2, tidyverse, lme4, MuMIn)
set.seed(123)

# Data Import
source("R/01_load-data_functions.R")
source("R/02_clean-data_functions.R")

# Linear Model
lm00 <- lm(Flower ~ Length + Width + Length_sq + Width_sq, data = inflo_data, na.action = "na.omit")
summary(lm00)
anova(lm00) # Length, Width, Width_sq is significant

# Calculating flower estimates based on LM
inflo_mixed <- combine_lm(inflo_data, lm00) 

# total number of flowers per individual
inflo_sum <- inflo_mixed %>%
  group_by(PlantID) %>%
  summarise(mixed = sum(mixed))


# Calculating ES indices
plant_data <- plant_yield(plant_traits, seed_data, inflo_sum) %>%
  replace_site()

# Export
write.csv(inflo_mixed, file = "data/processed/inflo_dat_LM.csv")
