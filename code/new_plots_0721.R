library(ggplot2)
library(tidyverse)
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

classifyStability <- Vectorize(function(stability){
  if (stability > 10){
    return("High")
  } else if (stability <= 10){
    return("Low")
  } else {
    return("Medium")
  }
})

dat <- read_csv('../results/060721_force03and01.csv') %>% 
  filter(collapsed!=TRUE) %>%
  # distinct(force, d, epsilon, stab, ratio, .keep_all=TRUE) %>%
  mutate(
    treatStab = classifyStability(stab),
    treatRat = classifyRatio(ratio)
  )

plotDat <- dat %>% 
  # filter(force==0|(force!=0&err<1E-6)) %>%
  group_by(d,epsilon, force, stab, ratio, fixed) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T),
    avgR = mean(avgR, na.rm=T)
  ) %>%
  filter(epsilon%in%c(1,5,10),
         ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15))

plotDat1 <- filter(plotDat, force!=0)

# Plot relatedness with respect to Y
ggplot(plotDat1) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gY,  color=as.factor(stab), linetype=as.factor(ratio)), size=1) +
  my_theme
# Plot relatedness with respect to X
ggplot(plotDat1) +
  facet_grid(~epsilon, scales="free") +
  geom_path(aes(avgR, gX,  color=as.factor(stab), linetype=as.factor(ratio)), size=1) +
  my_theme

#combined plot
plotDat <- dat %>% 
  filter(force!=0) %>%
  group_by(d,epsilon, force, stab, ratio, fixed) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T),
    avgR = mean(avgR, na.rm=T)
  ) %>%
  gather("trait", "val", gX, gY) %>%
  filter(epsilon==5,
         ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15))
ggplot(plotDat) +
  facet_grid(~trait, scales="free") +
  geom_path(aes(avgR, val,  color=as.factor(stab), linetype=as.factor(ratio)), size=1) +
  my_theme
ggsave("../graphs/traits_relatedness.pdf")

#combined plot removing stability 
plotDat <- dat %>% 
  filter(force!=0) %>%
  group_by(d,epsilon, force, ratio, fixed) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T),
    avgR = mean(avgR, na.rm=T)
  ) %>%
  gather("trait", "val", gX, gY) %>%
  filter(epsilon==5)
ggplot(plotDat) +
  facet_grid(~trait, scales="free") +
  geom_path(aes(avgR, val,  color=as.factor(ratio)), size=1) +
  my_theme
ggsave("../graphs/traits_allRatios.pdf")
# looking at stability and ratio
plotDat <- dat %>% 
  filter(force==0.03) %>%
  group_by(force, ratio,stab) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T)
  ) %>%
  gather("trait", "val", gX, gY) 
ggplot(plotDat) +
  facet_grid(~trait, scales="free") +
  geom_path(aes(ratio, val,  color=as.factor(stab)), size=1) +
  my_theme
ggsave("../graphs/traits_ratio.pdf")

plotDat <- dat %>% 
  filter(force==0.03) %>%
  group_by(force, stab, ratio) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T)
  ) %>%
  gather("trait", "val", gX, gY) 
ggplot(plotDat) +
  facet_grid(~trait, scales="free") +
  geom_path(aes(stab, val,  color=as.factor(ratio)), size=1) +
  my_theme
ggsave("../graphs/traits_stability.pdf")

# encounter rate
plotDat <- dat %>% 
  filter(force==0.03) %>%
  group_by(epsilon,force, ratio, d) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T)
  ) %>%
  gather("trait", "val", gX, gY)%>%
  filter(d==0.1)

ggplot(plotDat) +
  facet_grid(~trait, scales="free") +
  geom_path(aes(epsilon, val,  color=as.factor(ratio)), size=1) +
  my_theme

