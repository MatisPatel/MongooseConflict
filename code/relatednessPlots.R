library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
my_theme <- theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

dat <- read_csv('../results/third.csv') %>% 
  distinct(force, d, epsilon, gain, loss, .keep_all=TRUE) %>%
  mutate(gain_loss = paste(gain, loss, sep="_"),
         envRatio = gain/loss,
         avgRat = 1/((1/gain + 1/loss)/2),
         stability = 1/avgRat) 

classifyRatio <- Vectorize(function(envRatio){
  if (envRatio > 1){
    return("Bountiful")
  } else if (envRatio < 1){
    return("Desolate")
  } else {
    return("Neutral")
  }
})

classifyStability <- Vectorize(function(stability){
  if (stability > 10){
    return("High")
  } else if (stability < 10){
    return("Low")
  } else {
    return("Medium")
  }
})

plotDat <- dat %>% 
  filter(d>0, force!=1) %>%
  group_by(d, epsilon, force, envRatio, Stability) %>%
  summarise(
    gX = mean(groupAvgX),
    gY = mean(groupAvgY),
    avgR = mean(avgR),
    avgMort = mean(avgMort),
    avgFit = mean(meanFit),
    qVal = mean(qVal)
  )
