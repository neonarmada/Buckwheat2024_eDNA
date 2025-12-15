library(sf)
library(ggplot2)
library(tidyverse)
library(lme4)

inflo_data <- read.csv("data/inflo_dat.csv")
head(inflo_data)
inflo_data1 <- na.omit(inflo_data)

lm00<-lm(Flower~Length+Width+Length_sq+Width_sq,data=inflo_data1)
summary(lm00)
anova(lm00)
m00<-dredge(lm00,rank="AICc")
coef_00 <- get.models(m00, subset = 1)[[1]] %>%
  coefficients()
lm_data00 <- inflo_data %>%
  mutate(estimate00 = coef_00[[1]] + Length*coef_00[[2]] + Width*coef_00[[3]] +
           Width_sq*coef_00[[4]])

inflo_mixed <- lm_data00 %>%
  mutate(mixed = ceiling(coalesce(lm_data00$Flower,lm_data00$estimate00)))

inflo_sum <- inflo_mixed %>%
  group_by(PlantID) %>%
  summarise(flower = sum(mixed))

