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
  if (envRatio > 0.5){
    return("Bountiful")
  } else if (envRatio <= 0.5){
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
  } else if (stability <= 10){
    return("Low")
  } else {
    return("Medium")
  }
})


dat <- read_csv('../results/k02_run_full.csv') %>% 
  filter(ratio>0.2, fixed=="12", collapsed==FALSE) %>%
  distinct(fixed, force, d, epsilon, stab, ratio, .keep_all=TRUE) %>%
  mutate(gain_loss = paste(gain, loss, sep="_"),
         envRatio = gain/(gain+loss),
         stability = 1/(gain+loss)) %>%
  mutate(
    treatStab = classifyStability(stab),
    treatRat = classifyRatio(ratio)
  )


plotDat <- dat %>% 
  # filter(force==0|(force!=0&err<1E-6)) %>%
  group_by(d,epsilon, force, treatStab, ratio, fixed) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T),
    AgY = mean(AgroupAvgY, na.rm=T),
    SgY = mean(SgroupAvgY, na.rm=T),
    avgR = mean(avgR, na.rm=T),
    avgMort = mean(avgMort, na.rm=T),
    avgFit = mean(meanFit, na.rm=T),
    qVal = mean(qVal, na.rm=T),
    fit1 = mean(fit1, na.rm=T),
    fit2 = mean(fit2, na.rm=T),
    qVal1 = mean(diff(fit1, na.rm=T)),
    qVal2 = mean(diff(fit2, na.rm=T)),
    popSize = mean(popSize),
    err = mean(err)
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGX = gX/(gX+gY),
    # environment = paste(ratio, stab, sep="_")
  )%>%
  ungroup()

plotDat0 <- filter(plotDat, force==0)
plotDat5 <- filter(plotDat, force==0.1)

ggplot(filter(plotDat5, epsilon%in%c(1, 5, 10))) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gY,  color=treatStab, linetype=as.factor(ratio)), size=1) +
  my_theme
ggsave("../graphs/1_Y_rel_byEpsilon.png")

ggplot(filter(plotDat5, epsilon%in%c(1, 5, 10))) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gX,  color=treatStab, linetype=as.factor(ratio)), size=1) +
  my_theme

ggplot(filter(plotDat5, epsilon==10)) +
  geom_path(aes(d, gY,  color=as.factor(treatStab), linetype=as.factor(ratio)), size=1) +
  my_theme
ggplot(filter(plotDat0, epsilon==10)) +
  geom_path(aes(d, gY,  color=as.factor(stab), linetype=as.factor(ratio)), size=1) +
  my_theme

dd <- filter(plotDat5, 
             (ratio==0.2)|(ratio==0.1),
             epsilon==1
       )
ggplot(dd) +
  geom_path(aes(d, popSize,  color=as.factor(stab), linetype=as.factor(ratio)), size=1) +
  my_theme

ggplot(filter(plotDat5, epsilon%in%c(1,10), d>=0.1)) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gX,  color=treatStab, linetype=treatRat), size=1) +
  my_theme
ggsave("../graphs/2_X_rel_byEpsilon.png")

ggplot(filter(plotDat5, epsilon%in%c(1,10), d>0.1)) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gY/(gX+gY),  color=treatStab, linetype=treatRat), size=1) +
  my_theme
ggsave("../graphs/3_investY_rel_byEpsilon.png")

ggplot(filter(plotDat5,  epsilon%in%c(1,5,10), d>0.1)) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gX,  color=treatStab, linetype=as.factor(ratio)), size=1) +
  geom_path(aes(avgR, gY,  color=treatStab, linetype=as.factor(ratio)), size=1) +
  my_theme

ggplot(filter(plotDat5,  d%in%c(0.1, 0.5, 0.9))) +
  facet_grid(~d, scales="free") +
  geom_path(aes(epsilon, gX,  color=treatStab, linetype=treatRat), size=1) +
  geom_path(aes(epsilon, gY,  color=treatStab, linetype=treatRat), size=1) +
  my_theme

ggplot(filter(plotDat5,  d%in%c(0.1, 0.9))) +
  facet_grid(~d, scales="free") +
  geom_path(aes(epsilon, gY,  color=treatStab, linetype=treatRat), size=1) +
  my_theme
ggsave("../graphs/4_Y_epsilon_byD.png")

ggplot(filter(plotDat5,  d%in%c(0.5, 0.9))) +
  facet_grid(~d, scales="free") +
  geom_path(aes(epsilon, gX,  color=treatStab, linetype=treatRat), size=1) +
  my_theme
ggsave("../graphs/5_X_epsilon_byD.png")

ggplot(filter(plotDat5,  d%in%c(0.1, 0.9))) +
  facet_grid(~d, scales="free") +
  geom_path(aes(epsilon, gY/(gY+gX),  color=treatStab, linetype=treatRat), size=1) +
  my_theme
ggsave("../graphs/6_investY_epsilon_byD.png")

ggplot(filter(plotDat5,  d%in%c(0.1, 0.5, 0.9))) +
  facet_grid(~d, scales="free") +
  geom_path(aes(epsilon, popSize,  color=treatStab, linetype=treatRat), size=1) +
  my_theme

plotDat <- dat %>% 
  # filter(d>0, epsilon<20, force==0|(force==0.5&err<1E-7)) %>%
  group_by(force, stab, ratio) %>%
  summarise(
    gX = mean(groupAvgX),
    gY = mean(groupAvgY),
    AgY = mean(AgroupAvgY),
    SgY = mean(SgroupAvgY),
    avgR = mean(avgR),
    avgMort = mean(avgMort),
    avgFit = geomMean(meanFit),
    qVal = mean(qVal),
    fit1 = mean(fit1),
    fit2 = mean(fit2),
    qVal1 = mean(diff(fit1)),
    qVal2 = mean(diff(fit2))
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGX = gX/(gX+gY)
    # environment = paste(treatRat, treatStab, sep="_")
  )%>%
  ungroup()

plotDat0 <- filter(plotDat, force==0)
plotDat5 <- filter(plotDat, force==0.1)

dtemp <- filter(plotDat5) %>% 
  mutate(treatRat = classifyRatio(ratio)) %>% 
  group_by(stab, treatRat) %>%
  summarise(gY = mean(gY), gX = mean(gX))
ggplot(dtemp) +
  geom_path(aes((log(stab)), gY,  color=as.factor(treatRat)), size=1) +
  geom_path(aes((log(stab)), gX,  color=as.factor(treatRat)), size=1, linetype="dashed") +
  my_theme
ggsave("../graphs/7_Y_ratio.png")

dtemp <- filter(plotDat5) %>% 
  mutate(treatStab = classifyStability(stab)) %>% 
  group_by(ratio, treatStab) %>%
  summarise(gY = mean(gY), gX = mean(gX))
ggplot(dtemp) +
  geom_path(aes(ratio, gY,  color=as.factor(treatStab)), size=1) +
  geom_path(aes(ratio, gX,  color=as.factor(treatStab)), size=1, linetype="dashed") +
  my_theme
ggsave("../graphs/8_Y_stab.png")

dtemp <- filter(plotDat5) %>% 
  mutate(treatRat = classifyRatio(ratio)) %>% 
  group_by(stab, treatRat) %>%
  summarise(gX = mean(gX))
ggplot(dtemp) +
  geom_path(aes((log(stab)), gX,  color=as.factor(treatRat)), size=1) +
  geom_point(aes((log(stab)), gX,  color=as.factor(treatRat)), size=3) +
  my_theme
ggsave("../graphs/9_X_ratio.png")

dtemp <- filter(plotDat5) %>% 
  mutate(treatStab = classifyStability(stab)) %>% 
  group_by(ratio, treatStab) %>%
  summarise(gX = mean(gX))
ggplot(dtemp) +
  geom_path(aes(ratio, gX,  color=as.factor(treatStab)), size=1) +
  geom_point(aes(ratio, gX,  color=as.factor(treatStab)), size=3) +
  my_theme
ggsave("../graphs/10_X_stab.png")

ggplot(filter(plotDat5)) +
  geom_path(aes(stab, qVal,  color=as.factor(ratio)), size=1) +
  my_theme
ggplot(filter(plotDat5)) +
  geom_path(aes(stab, fit2,  color=as.factor(ratio)), size=1) +
  my_theme

ggplot(filter(plotDat5)) +
  geom_path(aes(stab, qVal,  color=as.factor(ratio)), size=1) +
  my_theme

ggplot(filter(plotDat5)) +
  geom_path(aes(ratio, qVal,  color=as.factor(stab)), size=1) +
  my_theme

ggplot(filter(plotDat0)) +
  geom_path(aes(stab, qVal,  color=as.factor(ratio)), size=1) +
  my_theme

ggplot(filter(plotDat0)) +
  geom_path(aes(ratio, qVal,  color=as.factor(stab)), size=1) +
  my_theme

ggplot(plotDat0) +
  geom_path(aes(ratio, qVal))
ggplot(plotDat5) +
  geom_path(aes(ratio, qVal))

ggplot(plotDat5) +
  geom_path(aes(ratio, gX)) +
  geom_path(aes(ratio, gY, linetype="Dashed"))
