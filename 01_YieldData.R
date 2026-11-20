p_load(sf, ggplot2, tidyverse, lme4, MuMIn)
set.seed(123)

# Data Import
inflo_data <- read.csv("data/raw/inflo_dat.csv")
head(inflo_data)

# Linear Model
lm00 <- lm(Flower ~ Length + Width + Length_sq + Width_sq, data = inflo_data, na.action = "na.omit")
summary(lm00)
anova(lm00)
coef_00 <- lm00 %>%
  coefficients()

# Adding flower estimates into
inflo_lm <- inflo_data %>%
  mutate(estimate00 = coef_00[[1]] + Length * coef_00[[2]] + Width * coef_00[[3]] +
    Width_sq * coef_00[[5]]) # Length_sq excluded bcs not significant

# note: some infloresence have negative amt of flowers
inflo_mixed <- inflo_lm %>%
  mutate(mixed = round(coalesce(lm_data00$Flower, lm_data00$estimate00)))

# total number of flowers per individual
inflo_sum <- inflo_mixed %>%
  group_by(PlantID) %>%
  summarise(flower = sum(mixed))

# Export
write.csv(inflo_mixed, file = "data/inflo_dat_LM.csv")
