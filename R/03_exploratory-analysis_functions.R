p_load(piecewiseSEM, reformulas)

# Convert to dummy variables: Variety, Type, Morph
dummy_VTM <- function(x) {
  Kitawase <- ifelse(x$Variety == "Kitawase", 1, 0)
  Leaf <- ifelse(x$Type == "Leaf", 1, 0)
  Pin <- ifelse(x$Morph == "Pin", 1, 0)
  
  x %>%
    mutate(Kitawase = Kitawase) %>%
    mutate(Leaf = Leaf) %>%
    mutate(Pin = Pin)
}

# eDNA model specification
eDNA_model <- function(x) {
  psem(
  lm(HerbivRate ~ NMDS1 + NMDS2, data = x),
  lm(SeedSet ~ NMDS1 + NMDS2, data = x),
  lm(DevRate ~ Variety + NMDS1 + NMDS2, data = x),
  lm(AvgMass ~ Variety + NMDS1 + NMDS2, data = x),
  lm(NMDS1 ~ ForestArea + Type + Variety + Morph, data = x),
  lm(NMDS2 ~ ForestArea + Type + Variety + Morph, data = x),
  x
  )
}

# eDNA_model <- '
# HerbivRate ~ Kitawase + NMDS1 + NMDS2
# SeedSet ~ Kitawase + NMDS1 + NMDS2
# DevRate ~ Kitawase + NMDS1 + NMDS2
# AvgMass ~ Kitawase + NMDS1 + NMDS2
# 
# NMDS1 ~ ForestArea + Leaf + Kitawase + Pin
# NMDS2 ~ ForestArea + Leaf + Kitawase + Pin
# '
