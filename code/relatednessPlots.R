library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
my_theme <- theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

classifyRatio <- Vectorize(function(envRatio){
  if (envRatio > 1){
    return("Bountiful")
  } else if (envRatio < 1){
    return("Desolate")
  } else {
    return("Neutral")
  }
})

harmMean <- function(vec){
  return(1/mean(1/vec))
}

geomMean <- function(vec){
  return(prod(vec)^(1/length(vec)))
}

classifyStability <- Vectorize(function(stability){
  if (stability > 10){
    return("High")
  } else if (stability < 10){
    return("Low")
  } else {
    return("Medium")
  }
})


dat <- read_csv('../results/second.csv') %>% 
  distinct(force, d, epsilon, gain, loss, .keep_all=TRUE) %>%
  mutate(gain_loss = paste(gain, loss, sep="_"),
         envRatio = gain/loss,
         avgRat = 1/((1/gain + 1/loss)/2),
         stability = 1/avgRat) %>%
  mutate(
    treatStab = classifyStability(stability),
    treatRat = classifyRatio(envRatio)
  )


plotDat <- dat %>% 
  filter(d>0, epsilon<20, force==0|(force==0.5&err<1E-7)) %>%
  group_by(d, epsilon, force, treatRat, treatStab) %>%
  summarise(
    gX = mean(groupAvgX),
    gY = mean(groupAvgY),
    AgY = mean(AgroupAvgY),
    SgY = mean(SgroupAvgY),
    avgR = mean(avgR),
    avgMort = mean(avgMort),
    avgFit = geomMean(meanFit),
    qVal = mean(qVal),
    fit1 = geomMean(fit1),
    fit2 = geomMean(fit2),
    qVal1 = mean(diff(fit1)),
    qVal2 = mean(diff(fit2))
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGX = gX/(gX+gY),
    environment = paste(treatRat, treatStab, sep="_")
  )%>%
  ungroup()

plotDat0 <- filter(plotDat, force==0)
plotDat5 <- filter(plotDat, force==0.5)

ggplot(filter(plotDat5, epsilon%in%c(1,5,10))) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, SgY,  color=environment), size=1) +
  my_theme
