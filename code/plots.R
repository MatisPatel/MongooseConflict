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
# %>%
#   filter(err<0.1)

plotDat <- filter(dat, 
                  d>0, epsilon<20) %>%
  group_by(tq, d, epsilon, treat) %>%
  summarise(
    gX = mean(groupAvgX),
    iX = mean(indAvgX),
    gY = mean(groupAvgY),
    iY = mean(indAvgY),
    AgX = mean(AgroupAvgX),
    AiX = mean(AindAvgX),
    AgY = mean(AgroupAvgY),
    AiY = mean(AindAvgY),
    SgX = mean(SgroupAvgX),
    SiX = mean(SindAvgX),
    SgY = mean(SgroupAvgY),
    SiY = mean(SindAvgY),
    avgR = mean(avgR),
    avgMort = mean(avgMort),
    fit1 = mean(fit1),
    fit2 = mean(fit2),
    avgFit = mean(meanFit)
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGYA = AgY/(AgX+AgY),
    investGYS = SgY/(SgX+SgY),
    investGX = gX/(gX+gY),
    investIY = iY/(iX+iY),
    totalG = gY+gX,
    totalGA = (AgX+AgY),
    totalGA = (SgX+SgY)
    )%>%
  ungroup()

#total Y
ggplot(filter(plotDat, epsilon %in% c(1,5,10))) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, investGY,  color=as.factor(treat)), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

ggplot(plotDat) +
  facet_grid(~treat, scales="free") +
  geom_point(aes((epsilon), gY/(mean(gY))), colour="green")+
  geom_point(aes((epsilon), gX/mean(gX))) +
  geom_line(aes(epsilon, gX/mean(gX)),
    stat="smooth",
    method ="lm", 
    se=FALSE,
    size =1) + 
  geom_line(aes(epsilon, gY/mean(gY)),
            stat="smooth",
            method ="lm", 
            se=FALSE,
            size =1, colour="green") +
  my_theme


# soical Y
ggplot(filter(plotDat)) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, investGYS,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

#asocial Y
ggplot(filter(plotDat)) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, investGYA,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

p2<- ggplot(plotDat) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, SgY,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

p1<- ggplot(plotDat) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, SgX,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme
grid.arrange(p1, p2, nrow = 1)
ggplot(plotDat) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, totalG,  color=treat), size=1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

fullplotDat <- filter(dat, 
                      d>0, epsilon<20) %>%
  group_by(d,epsilon, stability, envRatio) %>%
  summarise(
    gX = mean(groupAvgX),
    iX = mean(indAvgX),
    gY = mean(groupAvgY),
    iY = mean(indAvgY),
    AgX = mean(AgroupAvgX),
    AiX = mean(AindAvgX),
    AgY = mean(AgroupAvgY),
    AiY = mean(AindAvgY),
    SgX = mean(SgroupAvgX),
    SiX = mean(SindAvgX),
    SgY = mean(SgroupAvgY),
    SiY = mean(SindAvgY),
    avgR = mean(avgR),
    avgMort = mean(avgMort),
    fit1 = mean(fit1),
    fit2 = mean(fit2),
    avgFit = mean(meanFit)
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGYA = AgY/(AgX+AgY),
    investGYS = SgY/(SgX+SgY),
    investGX = gX/(gX+gY),
    investIY = iY/(iX+iY),
    totalG = gY+gX,
    totalGA = (AgX+AgY),
    totalGA = (SgX+SgY)
  )%>%
  ungroup()

pd2 <- fullplotDat %>%
  pivot_longer(c(investGYS, investGYA), names_to = "Y_ratio", values_to = "yval")

ggplot(pd2, aes(log(envRatio), avgFit)) +
  geom_point()

ggplot(pd2, aes(log(envRatio), yval)) +
  geom_point(aes(shape=Y_ratio, colour=as.factor(epsilon)), alpha=0.1) +
  geom_line(aes(linetype=Y_ratio, colour=as.factor(epsilon)),
            stat="smooth",
            method ="lm", 
            se=FALSE,
            size =1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

pd2 <- fullplotDat %>%
  pivot_longer(c(investGYS, investGYA), names_to = "Y_ratio", values_to = "yval")

ggplot(pd2, aes(log(envRatio), yval)) +
  geom_point(aes(shape=Y_ratio, colour=as.factor(epsilon)), alpha=0.1) +
  geom_line(aes(linetype=Y_ratio, colour=as.factor(epsilon)),
            stat="smooth",
            method ="lm", 
            se=FALSE,
            size =1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

ggplot(pd2, aes(log(stability), yval)) +
  geom_point(aes(shape=Y_ratio, colour=as.factor(epsilon)), alpha=0.1) +
  geom_line(aes(linetype=Y_ratio, colour=as.factor(epsilon)),
            stat="smooth",
            method ="lm", 
            se=FALSE,
            size =1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

fullplotDat <- filter(dat, 
                      d>0, epsilon<20) %>%
  group_by(epsilon, envRatio) %>%
  summarise(
    gX = mean(groupAvgX),
    iX = mean(indAvgX),
    gY = mean(groupAvgY),
    iY = mean(indAvgY),
    AgX = mean(AgroupAvgX),
    AiX = mean(AindAvgX),
    AgY = mean(AgroupAvgY),
    AiY = mean(AindAvgY),
    SgX = mean(SgroupAvgX),
    SiX = mean(SindAvgX),
    SgY = mean(SgroupAvgY),
    SiY = mean(SindAvgY),
    avgR = mean(avgR)
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGYA = AgY/(AgX+AgY),
    investGYS = SgY/(SgX+SgY),
    investGX = gX/(gX+gY),
    investIY = iY/(iX+iY),
    totalG = gY+gX,
    totalGA = (AgX+AgY),
    totalGA = (SgX+SgY)
  )%>%
  ungroup()

ggplot(fullplotDat) + 
  geom_path(aes((epsilon), gY/mean(gY), colour=as.factor(stability)))+
  geom_path(aes((epsilon), gX/mean(gX), colour=as.factor(stability)))

ggplot(fullplotDat) + 
  geom_path(aes((epsilon), gY, colour=as.factor(envRatio)))+
  geom_path(aes((epsilon), gX, colour=as.factor(envRatio)))

# ggplot(filter(fullplotDat)) +
#    geom_point(aes(log(stability), investGYS, colour=as.factor(epsilon)), 
#               alpha=0.1, shape=15) +
#      geom_point(aes(log(stability), investGYA, colour=as.factor(epsilon)), 
#               alpha=0.1, shape=17) + 
#      geom_line(aes(log(stability), investGYS, colour=as.factor(epsilon)),
#                stat="smooth",
#                method ="lm", 
#                se=FALSE,
#                size =1) +
#      geom_line(aes(log(stability), investGYA, colour=as.factor(epsilon)),
#                stat="smooth",
#                method ="lm", 
#                se=FALSE,
#                size =1, linetype="dashed") +
#      scale_colour_brewer(palette = "Dark2") +
#      my_theme

