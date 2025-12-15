library(sf)
library(ggplot2)
library(tidyverse)
library(lme4)
library(MuMIn)
library(ggsignif)
library(vegan)
set.seed(123)

##Importing Data & Prep
insect_raw <- read.csv("data/InsectSortingSheet_v5.csv")
head(insect_raw)
insect <- insect_raw %>%
  column_to_rownames(var="Site")

##NMDS
dist <- vegdist(insect, method = "bray")
nmds <- metaMDS(insect, autotransform=F)
plot(nmds)
nmds_result <- as.data.frame(scores(nmds,display=c("site")))
nmds_result

plot(nmds)
sp_scrs <- as.data.frame(scores(nmds,display=c("species")))
sp_scrs
taxon <- read.csv("data/insectsorting_taxon.csv")
sp_scrs2 <- sp_scrs %>%
  mutate(Order = taxon$Order)

##Envfit & Vectors
env <- read.csv("data/CollectedData_v4.csv")
env_site <- env %>%
  select(Site, Variety, ForestArea, SiteYield, HerbivRate) %>%
  distinct(Site, .keep_all = T) 
env1 <- env_site %>%
  slice(rep(1:n(), each = 2)) #duplicate each row to match trap replicates
env2 <- env1 %>%
  select(-Site, -SiteYield) 
vectors1 <- envfit(nmds,env2,perm=9999)
vectors1

##PERMANOVA
permanova.1 <- adonis2(insect~env2$ForestArea + env2$Variety,perm=9999,by="terms")
#write.csv(permanova.1, "output/TrapNMDS_permanova_envfactor.csv")

centroids <- scores(vectors, display = "factors") %>%
  data.frame() %>%
  rownames_to_column(var = "Variety")
centroids <- centroids %>%
  mutate(Variety = str_remove(centroids$Variety, "Variety"))

vectors_scores <- as.data.frame(scores(vectors, display = "vectors"))
vectors_scores$factors <- rownames(vectors_scores)
vectors_scores

## Can trap NMDS explain changes in yield and herbivory rate?
nmds_combined <- nmds_result %>%
  bind_cols(env1)

lm_yield <- lm(SiteYield ~ NMDS1 + NMDS2, nmds_combined)
summary(lm_yield)
anova(lm_yield)

lm_herb <- lm(HerbivRate ~ NMDS1 + NMDS2, nmds_combined)
summary(lm_herb)
anova(lm_herb)

# Graphics ======
# Site Plot
nmds_site <- ggplot(nmds_result, aes(x = NMDS1, y = NMDS2)) +
  #env_vectors
  geom_segment(data = vectors_scores, aes(x=0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  #species scores (square)
  #geom_point(data=sp_scrs2, aes(fill = Order), size = 3, shape = 22, stroke= 1)+
  #Site Data (Circle)
  geom_point(shape =21, size = 8, stroke = 1, aes(fill = env$Site)) + 
  #Centroids (triangle)
  geom_point(data = centroids, shape = 24, aes(fill = Variety), size = 6, stroke = 1) + 
  #colours
  scale_fill_manual(values = c("A"="white", "B"="#333333","C"="#666666", 
                               "D"= "#CCCCCC", "E" = "#999999", "F" = "black",
                               "Horominori" = "#9A287A", "Kitawase" = "#CC949D",
                               "Diptera" = "#ab384a", "Hymenoptera" = "#ffe699",
                               "Other" = "#91bcae"))+
  geom_text(data=vectors_scores, aes(x=NMDS1,y=NMDS2, label=factors),size=6) +
  theme_bw(base_size = 28) +
  theme(legend.position = "none", axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
nmds_site
ggsave("output/trapnmds_site.png", 
       plot = nmds_site, width = 2400 , height = 1600, units = "px", dpi = 300)

#species plot
nmds_sp <- ggplot(nmds_result, aes(x = NMDS1, y = NMDS2)) +
  #env_vectors
  geom_segment(data = vectors_scores, aes(x=0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  #species scores (square)
  geom_point(data=sp_scrs2, aes(fill = Order), size = 3, shape = 22, stroke= 1)+
  #Centroids (triangle)
  geom_point(data = centroids, aes(fill = Variety), shape = 24, size = 6, stroke = 1)+
  #colours
  scale_fill_manual(values = c("A"="white", "B"="#333333","C"="#666666", 
                               "D"= "#CCCCCC", "E" = "#999999", "F" = "black",
                               "Horominori" = "#9A287A", "Kitawase" = "#CC949D",
                               "Diptera" = "#ab384a", "Hymenoptera" = "#ffe699",
                               "Other" = "#91bcae"))+
  geom_text(data=vectors_scores, aes(x=NMDS1,y=NMDS2, label=factors),size=6) +
  theme_bw(base_size = 28) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"), legend.position = "none",
        legend.justification = c("right", "top"), legend.text = element_text(size = 18),
        legend.background = element_rect(colour = "black"))
nmds_sp
ggsave("output/trapnmds_sp.png", 
       plot = nmds_sp, width = 2400 , height = 1600, units = "px", dpi = 300)

