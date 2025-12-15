library(sf)
library(ggplot2)
library(tidyverse)
library(lme4)

inflo_data <- read.csv("data/inflo_data.csv")
head(inflo_data)
inflo_data1 <- na.omit(inflo_data)

lm00<-lm(Flower~Length+Width+Length_sq+Width_sq,data=inflo_data1)
m00<-dredge(lm00,rank="AICc")
coef_00 <- get.models(m00, subset = 1)[[1]] %>%
  coefficients()
lm_data00 <- inflo_data %>%
  mutate(estimate00 = coef_00[[1]] + Length*coef_00[[2]] + Width*coef_00[[3]] +
           Width_sq*coef_00[[4]])

#add all lm00-lm03 estimations to the same table
#compare via scatters

lm01<-lm(Flower~Length_sq+Width,data=inflo_data1)
summary(lm01)
anova(lm01)
coef_01 <- coefficients(lm01)
lm_data01 <- lm_data00 %>%
  mutate(estimate01 = coef_01[[1]] + Length_sq*coef_01[[2]] + Width*coef_01[[3]])

lm02<-lm(Flower~Length+Width_sq,data=inflo_data1)
summary(lm02)
anova(lm02)
coef_02 <- coefficients(lm02)
lm_data02 <- lm_data01 %>%
  mutate(estimate02 = coef_02[[1]] + Length*coef_02[[2]] + Width_sq*coef_02[[3]])

lm03<-lm(Flower~Length_sq+Width_sq,data=inflo_data1)
summary(lm03)
anova(lm03)
coef_03 <- coefficients(lm03)
lm_data03 <- lm_data02 %>%
  mutate(estimate03 = coef_03[[1]] + Length_sq*coef_03[[2]] + Width_sq*coef_03[[3]])

lm_compare <- lm_data03 %>%
  na.omit(Flower)

ggplot(data = lm_compare, aes(x = Flower)) +
  geom_point(aes(y = estimate00), colour = "red", alpha = 1) +
  geom_point(aes(y = estimate01), colour = "green", alpha = 0.1) +
  geom_point(aes(y = estimate02), colour = "blue", alpha = 0.1) +
  geom_point(aes(y = estimate03), colour = "orange", alpha = 0.1) +
  geom_abline(aes(intercept = 0, slope = 1), size = 1) +
  ylim(0, 500) +
  theme_classic()
