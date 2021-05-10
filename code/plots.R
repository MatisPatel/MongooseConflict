library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
my_theme <- theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

classifier <- function(gain, loss){
  if (gain == loss){
    return("neutral")
  } else if (gain > loss){
    return("benign")
  } else {
    return("harsh")
  }
}

vclassifier <- Vectorize(classifier)

dat <- read_csv('../results/firstRun') %>% 
  distinct(d, epsilon, gain, loss, .keep_all=TRUE) %>%
  mutate(gain_loss = paste(gain, loss, sep="_"),
         envRatio = gain/loss,
         avgRat = 1/((1/gain + 1/loss)/2),
         stability = 1/avgRat,
         treat = vclassifier(gain, loss))

plotDat <- filter(dat, 
                  tn==2,
                  d>0, epsilon<20) %>%
  group_by(d, epsilon, treat) %>%
  summarise(
    gX = mean(groupAvgX),
    iX = mean(indAvgX),
    gY = mean(groupAvgY),
    iY = mean(indAvgY),
    avgR = mean(R)
  ) %>%
  ungroup()

ggplot(plotDat) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(d, gY,  color=treat)) +
  my_theme

ggplot(plotDat) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(d, gY,  color=treat)) +
  my_theme

