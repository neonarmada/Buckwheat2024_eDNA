library(sf)
library(ggplot2)
library(tidyverse)
library(lme4)
library(MuMIn)
library(ggsignif)

buckwheat <- read.csv("data/CollectedData_v4.csv")
head(buckwheat)

## Mature Seed Set x Forest Area Plot
devrate_plot <- ggplot(data = buckwheat, aes(x = ForestArea, y = DevRate)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#E76D11") +
  labs( x = "% Forest Area", y = "Mature Seed Set") +
  theme_classic(base_size=28)
devrate_plot

## Total Seed Set x Forest
seedset_plot <- ggplot(data = buckwheat, aes(x = ForestArea, y = SeedSet)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#3F73BF") +
  labs( x = "% Forest Area", y = "Total Seed Set") +
  theme_classic(base_size=28)
seedset_plot

## Mature Seed Set x Variety Averages
dev_avg <- data.frame(
  Var = c("Horominori", "Kitawase"),
  Avg = c(0.0276, 0.0366666666666667))
## Variety x DevRate Plot
dev_var <- ggplot(data = buckwheat, aes(x = Variety, y = DevRate)) +
  geom_col(data = dev_avg, aes(x = Var, y = Avg, fill = Var), width = 0.3,
           colour = "black", linewidth = .6) +
  scale_fill_manual(values = c("Horominori" = "#9B287B", "Kitawase" = "#CC949D")) +
  geom_point() +
  ## Remove gap between axis and plot
  scale_y_continuous(expand = c(0,0)) +
  geom_signif(comparisons = list(c("Horominori","Kitawase")), annotations = "*", 
              y_position = 0.085, textsize = 10, vjust = 0.5) +
  coord_cartesian(ylim = c(0,0.095))+
  labs(y = "Mature Seed Set") +
  theme_classic(base_size=28) +
  theme(legend.position = "none", axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
dev_var

## Herb_Variety Averages For Bar Plot
herb_av <- data.frame(
  Var = c("Horominori", "Kitawase"),
  Avg = c(83.374, 59.85778))
## Variety x Herbivory Plot
herbiv_graph <- ggplot(data = buckwheat, aes(x = Variety, y = HerbivRate)) +
  geom_col(data = herb_av, aes(x = Var, y = Avg, fill = Var), width = 0.3,
           colour = "black", linewidth = 1) +
  scale_fill_manual(values = c("Horominori" = "#9B287B", "Kitawase" = "#CC949D")) +
  geom_point() +
  coord_cartesian(ylim = c(0,115)) +
  ## Remove gap between axis and plot
  scale_y_continuous(expand = c(0,0)) +
  geom_signif(comparisons = list(c("Horominori","Kitawase")), annotations = "***",
              textsize = 10, vjust = 0.5, y_position = 105) +
  labs(y = "Leaf Herbivory Rate") +
  theme_classic(base_size=28) +
  theme(legend.position = "none", axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
herbiv_graph

## Site Yield x Forest Area Plot
siteyield_plot <- ggplot(data = buckwheat, aes(x = ForestArea, y = SiteYield)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE)
siteyield_plot

## Site Yield vs Seed Set
ggplot(data = buckwheat, aes(x = ForestArea)) +
  geom_col(aes(y = SiteYield)) +
  geom_point(aes(y = DevRate * 2000)) +
  scale_y_continuous(name = "SiteYield (Bar)", 
                     sec.axis = sec_axis(transform = ~./2000, 
                                         name = "DevRate (Scatter)")) +
  theme_classic()

## Abundance 
abun <- read.csv("data/InsectSortingSheet_abun.csv")
abun_arr <- arrange(abun, desc(Abun_r)) %>%
  mutate(Taxon = factor(Taxon, levels = Taxon))
head(abun_arr)
rel_abun <- ggplot(data = abun_arr, aes(x = Taxon, y = Abun_r, fill = Order)) +
  geom_col(color = "black", position = "dodge", width = 1) +
  scale_fill_manual(values = c("Diptera" = "#ab384a", "Hymenoptera" = "#91bcae", 
                               "Other" = "#b4c9e7")) +
  labs(y = "Rel. Abundance") +
  ##remove gaps between plot and axis
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  coord_cartesian(ylim = c(0,0.18))+
  theme_classic(base_size = 28) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, colour = "black"), 
        axis.text.y = element_text(colour = "black"), legend.position = "inside",
        legend.justification = c("right", "top"), legend.text = element_text(size = 18),
        legend.background = element_rect(colour = "black"))
rel_abun

##VarPar Plot, Relative importance
vartpart <- read.csv("output/VarPar_table.csv")
head(vartpart)
vartpart_arr <- arrange(vartpart, Factor)
#Arrange Factors so that residue is on the bottom
vartpart_arr$Factor <- factor(
  vartpart_arr$Factor, levels = c("ForestArea","Variety","Morph","HerbivRate","Residue"))

head(vartpart_arr)
#Stacked Elements plot (GreenPurplOrange)
vartpart_plot <- ggplot(data = vartpart_arr, aes(x= VarPar, y= percent, fill = Factor)) +
  geom_col(position = "fill", width = 0.25, colour = "black", linewidth = .55) +
  scale_fill_manual(values = c("Residue"= "darkgrey","Variety"="#9B287B","Morph" = "#63A8AF",
                               "HerbivRate"= "#AEC6C4", "ForestArea"="#FFC000"))+
  theme_void() +
  theme(legend.position = "none")
vartpart_plot

#VarPar, v2
vartpart3 <- read.csv("output/VarPar_table4.csv")
vartpart_arr3 <- arrange(vartpart3, Factor)
vartpart_arr3$Factor <- factor(
  vartpart_arr3$Factor, levels = c("ForestArea","Variety","Morph", 
                                   "ForestAreaMorph", "ForestAreaVariety", 
                                   "MorphVariety", "ForestAreaMorphVariety","Residuals"))
view(vartpart_arr3)
#Stacked Elements plot (GreenPurplOrange)
vartpart_plot <- ggplot(data = vartpart_arr3, aes(x= VarPar, y= percent, fill = Factor)) +
  geom_col(position = "fill", width = 0.3, colour = "black", linewidth = .75) +
  scale_fill_manual(name = "Relative Contribution", breaks =c("ForestArea", "Variety", 
                                                              "Morph", "GxE Interactions", "Residuals"),
                    values = c("Residuals"= "grey","Variety"="#9B287B","Morph" = "#63A8AF",
                               "ForestArea"="#FFC000", "GxE Interactions"="#7F7F7F"), )+
  theme_void() +
  theme(legend.position = "none")
vartpart_plot


# ggsave
ggsave("output/devrate_plot.png", plot = devrate_plot, 
       width = 2400 , height = 1600, units = "px", dpi = 300)
ggsave("output/seedset_plot.png", plot = seedset_plot, 
       width = 2400 , height = 1600, units = "px", dpi = 300)
ggsave("output/devvar_plot.png", plot = dev_var, 
       width = 2400 , height = 1600, units = "px", dpi = 300)
ggsave("output/herbiv_plot.png", plot = herbiv_graph, 
       width = 2400 , height = 1600, units = "px", dpi = 300)
ggsave("output/relabun_plot.png", plot = rel_abun, 
       width = 2400 , height = 1600, units = "px", dpi = 300)
ggsave("output/vartpart_plot.png", plot = vartpart_plot, 
       width = 1200 , height = 800, units = "px", dpi = 300)
ggsave("output/vartpart_plot2.png", plot = vartpart_plot, 
       width = 1500, height = 1100, units = "px", dpi = 300)