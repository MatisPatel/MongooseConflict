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

dat <- read_csv('../results/third.csv') %>% 
  distinct(force, d, epsilon, gain, loss, .keep_all=TRUE) %>%
  mutate(gain_loss = paste(gain, loss, sep="_"),
         envRatio = gain/loss,
         avgRat = 1/((1/gain + 1/loss)/2),
         stability = 1/avgRat,
         treat = vclassifier(gain, loss)) 
# %>%
#   filter(err<0.1)

dat2 <- dat %>%
  gather("patchW", "W", starts_with("tW")) %>%
  extract(patchW, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")

ggplot(filter(dat2, epsilon%in%c(1,5,10))) +
  facet_grid(vars(d), vars(epsilon)) +
  geom_bar(aes(tn, W, fill=tq), position = position_dodge(), stat = "identity") +
  theme_bw()

dat2 <- dat %>%
  gather("patchM", "M", starts_with(c("tY", "tX"))) %>%
  extract(patchM, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  spread(measure, M) %>%
  mutate(tM =  mortfun(tX, tY, as.numeric(tn)))

mortfun <- function(tX, tY, tn){
  return(0.1*2.718^(-1 * (tn)*((tX*(tn-1) + tX)/tn)) + 0.01*tX^2 + 0.01*tY^2)
}

ggplot(filter(dat2, epsilon%in%c(1,5,10), d%in%c(0.1, 0.5, 0.9), force==0.5)) +
  facet_grid(vars(d), vars(epsilon)) +
  geom_bar(aes(tn, tM, fill=tq), position = position_dodge(), stat = "identity") +
  theme_bw()

plotDat <- filter(dat, 
                  d>0, epsilon<20) %>%
  group_by(d, epsilon, treat, force) %>%
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
    avgFit = mean(meanFit),
    qVal = mean(qVal)
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
ggplot(filter(plotDat, epsilon %in% c(1,5,10), force==0.5)) +
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
  geom_path(aes(avgR, investGY,  color=treat), size=1) +
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
  group_by(force, d,epsilon, stability, envRatio) %>%
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
    avgFit = mean(meanFit),
    avgQ = mean(qVal)
  ) %>%
  mutate(
    investGY = gY/(gX+gY),
    investGYA = AgY/(AgX+AgY),
    investGYS = SgY/(SgX+SgY),
    investGX = gX/(gX+gY),
    investGXA = AgX/(AgX+AgY),
    investGXS = SgX/(SgX+SgY),
    investIY = iY/(iX+iY),
    totalG = gY+gX,
    totalGA = (AgX+AgY),
    totalGA = (SgX+SgY)
  )%>%
  ungroup()

pd2 <- fullplotDat %>%
  pivot_longer(c(investGYS, investGYA), names_to = "Y_ratio", values_to = "yval")

ggplot(filter(pd2, force!=1), aes(log(envRatio), diff)) +
  geom_point(aes(colour=force)) +
  geom_line(aes(log(envRatio), diff, linetype=as.factor(force)),
            stat="smooth",
            method ="lm", 
            se=FALSE,
            size =1)

ggplot(filter(pd2, force==0.5)) +
  geom_boxplot(aes(as.factor(log(stability)), 
                 fit1, 
                 colour="red")) +
  geom_boxplot(aes(as.factor(log(stability)), 
                 fit2, 
                 colour="black")) 

ggplot(filter(pd2, force==0.5), aes(log(envRatio), yval)) +
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

ggplot(filter(fullplotDat, force==0.5), aes(log(envRatio), log(avgQ+1))) +
  geom_point(aes(colour=as.factor(epsilon)), alpha=0.1) +
  geom_line(aes(colour=as.factor(epsilon)),
            stat="smooth",
            method ="lm", 
            se=FALSE,
            size =1) +
  scale_colour_brewer(palette = "Dark2") +
  my_theme

ggplot(filter(fullplotDat, force==0.5), aes(log(stability), log(avgQ+1))) +
  geom_point(aes(colour=as.factor(epsilon)), alpha=0.1) +
  geom_line(aes(colour=as.factor(epsilon)),
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

pd2 <- fullplotDat %>%
  pivot_longer(c(gX, gY), names_to = "Y_ratio", values_to = "yval")

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
                      d>0, epsilon<20, err<1E-9) %>%
  group_by(d, epsilon, stability) %>%
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
    avgFit = mean(meanFit),
    avgQ = mean(qVal)
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

ggplot(filter(fullplotDat, d==0.5)) + 
  geom_path(aes((epsilon), gY/mean(gY), colour=as.factor(stability)))+
  geom_path(aes((epsilon), gX/mean(gX), colour=as.factor(stability)))

ggplot(filter(fullplotDat, d==0.5)) + 
  geom_path(aes((epsilon), avgMort, colour=as.factor(stability)))+
  geom_path(aes((epsilon), avgMort, colour=as.factor(stability)))

# point corr plots 
ggplot(filter(fullplotDat, d==0.5)) +
  geom_path(aes(epsilon, avgQ, colour=as.factor(stability)))

fullplotDat <- filter(dat, 
                      d>0, epsilon<20, err<1E-9) %>%
  group_by(d, epsilon, envRatio) %>%
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
    avgFit = mean(meanFit),
    avgQ = mean(qVal)
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

ggplot(filter(fullplotDat, d==0.5)) + 
  geom_path(aes((epsilon), gY/mean(gY), colour=as.factor(envRatio)))+
  geom_path(aes((epsilon), gX/mean(gX), colour=as.factor(envRatio)))

ggplot(filter(fullplotDat, d==0.5)) + 
  geom_path(aes((epsilon), diff, colour=as.factor(envRatio)))+
  geom_path(aes((epsilon), diff, colour=as.factor(envRatio)))

# point corr plots 
ggplot(filter(fullplotDat, d==0.5)) +
  geom_path(aes(epsilon, avgQ, colour=as.factor(envRatio)))

ggplot(filter(fullplotDat, d==0.5)) + 
  geom_point(aes((log(envRatio)), gY/mean(gY), colour=as.factor(epsilon)))


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

