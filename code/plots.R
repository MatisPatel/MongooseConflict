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

dat <- read_csv('../results/firstRun.csv') %>% 
  distinct(d, epsilon, gain, loss, .keep_all=TRUE) %>%
  mutate(gain_loss = paste(gain, loss, sep="_"),
         envRatio = gain/loss,
         avgRat = 1/((1/gain + 1/loss)/2),
         stability = 1/avgRat,
         treat = vclassifier(gain, loss))

plotDat <- filter(dat, 
                  d>0, epsilon<20) %>%
  group_by(tn, d, epsilon, treat) %>%
  summarise(
    gX = mean(groupAvgX),
    iX = mean(indAvgX),
    gY = mean(groupAvgY),
    iY = mean(indAvgY),
    avgR = mean(avgR)
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGX = gX/(gX+gY),
    investIY = iY/(iX+iY),
    totalG = gY+gX
    )%>%
  ungroup()

# soical Y
ggplot(filter(plotDat, tn==2)) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, investGY,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

#asocial Y
ggplot(filter(plotDat, tn==1)) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, investGY,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

ggplot(plotDat) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gY,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

ggplot(plotDat) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gX,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

ggplot(plotDat) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, totalG,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

fullplotDat <- filter(dat, 
                      tn==2,
                      d>0, epsilon<20) %>%
  group_by(d,epsilon, stability, envRatio) %>%
  summarise(
    gX = mean(groupAvgX),
    iX = mean(indAvgX),
    gY = mean(groupAvgY),
    iY = mean(indAvgY),
    avgR = mean(avgR)
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGX = gX/(gX+gY),
    investIY = iY/(iX+iY),
    totalG = gY+gX
  )%>%
  ungroup()


ggplot(filter(fullplotDat), 
       aes(log(envRatio), investGY, colour=as.factor(epsilon))) +
  geom_point() + 
  geom_line(stat="smooth",
            method ="lm", 
            se=FALSE,
            size =1, alpha=0.4) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

ggplot(filter(fullplotDat), 
       aes(log(stability), investGY, colour=as.factor(epsilon))) +
  geom_point() + 
  geom_line(stat="smooth",
            method ="lm", 
            se=FALSE,
            size =1, alpha=0.4)+
  scale_colour_brewer(palette = "Dark2") +
  my_theme

