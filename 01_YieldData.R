library(sf)
library(ggplot2)
library(tidyverse)
library(lme4)
library(MuMIn)

inflo_data <- read.csv("data/inflo_dat.csv")
head(inflo_data)

lm00<-lm(Flower~Length+Width+Length_sq+Width_sq,data=inflo_data, na.action = "na.omit")
summary(lm00)
anova(lm00)
coef_00 <- lm00 %>%
  coefficients()
lm_data00 <- inflo_data %>%
  mutate(estimate00 = coef_00[[1]] + Length*coef_00[[2]] + Width*coef_00[[3]] +
           Width_sq*coef_00[[5]]) #Length_sq excluded bcs not significant

inflo_mixed <- lm_data00 %>%
  mutate(mixed = round(coalesce(lm_data00$Flower,lm_data00$estimate00)))

inflo_sum <- inflo_mixed %>%
  group_by(PlantID) %>%
  summarise(flower = sum(mixed))

#Export
write.csv(inflo_mixed, file = "data/inflo_dat_LM.csv")
