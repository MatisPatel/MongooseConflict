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
  group_by(force, ratio) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T)
  ) %>%
  gather("trait", "val", gX, gY) 
ggplot(plotDat) +
  facet_grid(~trait, scales="free") +
  geom_path(aes(ratio, val), size=1) +
  my_theme
ggsave("../graphs/traits_ratio.png")

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
ggsave("../graphs/traits_encounterRate.pdf")

plotDat <- dat %>% 
  filter(force!=0) %>%
  group_by(d, force, ratio, stab, fixed) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T),
    avgR = mean(avgR, na.rm=T),
    avgM = mean(avgMort, na.rm=T)
  ) %>%
  gather("trait", "val", gX, gY) %>%
  filter(
         force==0.03,
         ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15))
ggplot(plotDat) +
  facet_grid(~trait, scales="free") +
  geom_path(aes(avgR, val,  color=as.factor(ratio), linetype=as.factor(stab)), size=1) +
  my_theme
ggsave("../graphs/traits_relatedness.pdf")

plotDat <- dat %>% 
  filter(force!=0) %>%
  group_by(d, force, ratio, stab, fixed) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T),
    avgR = mean(avgR, na.rm=T),
    avgM = mean(avgMort, na.rm=T)
  ) %>%
  gather("trait", "val", gX, gY) %>%
  filter(
    force==0.03,
    ratio%in%c(0.1, 0.5, 0.9),
    stab%in%c(5, 10, 15))
ggplot(plotDat) +
  facet_grid(ratio~trait, scales="free") +
  geom_path(aes(val, avgM,  color=as.factor(ratio), linetype=as.factor(stab)), size=1) +
  my_theme
ggsave("../graphs/traits_mortality.pdf")

plotDat <- dat %>% 
  filter(force!=0) %>%
  group_by(d, force, ratio, stab, fixed) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T),
    avgR = mean(avgR, na.rm=T),
    avgM = mean(avgMort, na.rm=T)
  ) %>%
  # gather("trait", "val", gX, gY) %>%
  filter(
    force==0.03,
    ratio%in%c(0.1, 0.5, 0.9),
    stab%in%c(5, 10, 15))
ggplot(plotDat) +
  # facet_grid(~ratio, scales="free") +
  geom_path(aes(avgR, avgM,  color=as.factor(ratio), linetype=as.factor(stab)), size=1) +
  my_theme
ggsave("../graphs/relatedness_mortality.pdf")

plotDat <- dat %>% 
  filter(force!=0) %>%
  group_by(d, force, ratio, stab, fixed) %>%
  summarise(
    gX = mean(groupAvgX, na.rm=T),
    gY = mean(groupAvgY, na.rm=T),
    avgR = mean(avgR, na.rm=T),
    avgM = mean(avgMort, na.rm=T)
  ) %>%
  # gather("trait", "val", gX, gY) %>%
  filter(
    force==0.03,
    d%in%c(0.1, 0.5, 0.9),
    stab%in%c(5, 10, 15))
ggplot(plotDat) +
  # facet_grid(~ratio, scales="free") +
  geom_path(aes(ratio, avgM, colour=as.factor(stab),linetype=as.factor(d)), size=1) +
  my_theme
ggsave("../graphs/relatedness_mortality.pdf")

plotDat <- dat %>%
  gather("measure", "val", starts_with(c("tX", "tY"))) %>%
  select(measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, force, ratio, fixed, measure, tq, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T)
  ) %>%
  mutate(tq = as.numeric(tq),
         tn = as.numeric(tn)) %>%
  filter(measure%in%c("tX","tY"), 
         force==0.03,
         ratio%in%c(0.1, 0.5, 0.9),
         # stab%in%c(5, 10, 15),
         tn==2
         ) %>%
  ungroup()
  
ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(~measure, scales="free") +
  geom_path(aes(tq, meanVal, 
                colour=as.factor(ratio)),
            size=2) +
  my_theme
ggsave("../graphs/quality_trait.png")

ggplot(plotDat%>%filter(stab==10, 
                        d==0.5, 
                        epsilon==5,
                        tn==1)) + 
  geom_tile(aes(tq, ratio, fill=meanVal)) + 
  scale_fill_viridis_c()
  

plotDat <- dat %>% filter(force==0.03) %>%
  gather("measure", "val", starts_with(c("tR"), ignore.case = FALSE)) %>%
  select(avgMort, avgR, measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, ratio, stab, fixed, measure, tq, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T),
  ) %>%
  mutate(tq = as.numeric(tq),
         tn = as.numeric(tn)) %>%
  filter(force==0.03,
         ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         tn>0
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(~measure, scales="free") +
  geom_path(aes(tq, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/quality_R.pdf")  

plotDat <- dat %>% filter(force==0.03) %>%
  gather("measure", "val", starts_with(c("M"), ignore.case = FALSE)) %>%
  select(avgMort, avgR, measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, ratio, stab, fixed, measure, tq, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T),
  ) %>%
  mutate(tq = as.numeric(tq), tn = as.numeric(tn)) %>%
  filter(ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         measure%in%c("Mf"),
         tn>0
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(tn~measure, scales="free") +
  geom_path(aes(tq, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/quality_Mort.pdf")  

plotDat <- dat %>% filter(force==0.03) %>%
  gather("measure", "val", starts_with(c("relW"), ignore.case = FALSE)) %>%
  select(avgMort, avgR, measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{4})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, ratio, stab, fixed, measure, tq, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T),
  ) %>%
  mutate(tq = as.numeric(tq), tn=as.numeric(tn)) %>%
  filter(ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         tn>0
         # measure%in%c("Mf")
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(tn~measure, scales="free") +
  geom_path(aes(tq, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/quality_fitness.pdf")  

plotDat <- dat %>% filter(force==0.03) %>%
  gather("measure", "val", starts_with(c("M"), ignore.case = FALSE)) %>%
  select(avgMort, avgR, measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, ratio, stab, fixed, measure, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T),
  ) %>%
  mutate(tq = as.numeric(tn)) %>%
  filter(ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         measure%in%c("Mf"),
         tq>0
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(~measure, scales="free") +
  geom_path(aes(tq, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/size_Mort.pdf") 

plotDat <- dat %>% filter(force==0.03) %>%
  gather("measure", "val", starts_with(c("relW"), ignore.case = FALSE)) %>%
  select(avgMort, avgR, measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{4})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, ratio, stab, fixed, measure, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T),
  ) %>%
  mutate(tq = as.numeric(tn)) %>%
  filter(ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         tq>0
         # measure%in%c("Mf")
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(~measure, scales="free") +
  geom_path(aes(tq, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/size_fitness.pdf") 

### Data for 3x5
dat2 <- read_csv('../results/150721_53_and_35.csv') %>% 
  filter(collapsed!=TRUE, n==5, q==3)

plotDat <- dat2 %>% filter(force==0.03) %>%
  gather("measure", "val", starts_with(c("M"), ignore.case = FALSE)) %>%
  select(avgMort, avgR, measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, ratio, stab, fixed, measure, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T),
  ) %>%
  mutate(tn = as.numeric(tn)) %>%
  filter(ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         measure%in%c("Mf"),
         tn>0
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(~measure, scales="free") +
  geom_path(aes(tn, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme

ggsave("../graphs/size_Mort2.pdf") 

plotDat <- dat2 %>% filter(force==0.03) %>%
  gather("measure", "val", starts_with(c("relW"), ignore.case = FALSE)) %>%
  select(avgMort, avgR, measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{4})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, ratio, stab, fixed, measure, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T),
  ) %>%
  mutate(tn = as.numeric(tn)) %>%
  filter(ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         tn>0
         # measure%in%c("Mf")
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(~measure, scales="free") +
  geom_path(aes(tn, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/size_fitness2.pdf") 

plotDat <- dat2 %>%
  gather("measure", "val", starts_with(c("tX", "tY"))) %>%
  select(measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, force, ratio, stab, fixed, measure, tq, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T)
  ) %>%
  mutate(tq = as.numeric(tq),
         tn = as.numeric(tn)) %>%
  filter(measure%in%c("tX","tY"), 
         force==0.03,
         ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         tn>0, tq<=3
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(tq~measure, scales="free") +
  geom_path(aes(tn, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/size_trait2.pdf")

plotDat <- dat2 %>%
  gather("measure", "val", starts_with(c("tX", "tY"))) %>%
  select(measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  group_by(d, epsilon, force, ratio, stab, fixed, measure, tq, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T)
  ) %>%
  mutate(tq = as.numeric(tq),
         tn = as.numeric(tn)) %>%
  filter(measure%in%c("tX","tY"), 
         force==0.03,
         ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         tn>0, tq<=3
  ) %>%
  ungroup()

ggplot(plotDat %>% filter(d==0.5, epsilon==5)) + 
  facet_grid(tq~measure, scales="free") +
  geom_path(aes(tn, meanVal*tn, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/group_size_trait2.pdf")

plotDat <- dat2 %>% filter(force==0.03) %>%
  gather("measure", "val", starts_with(c("tR"), ignore.case = FALSE)) %>%
  select(avgMort, avgR, measure, val, d, epsilon, stab, ratio, force, fixed) %>%
  extract(measure, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)")%>%
  group_by(d,ratio, stab, fixed, measure, tn) %>%
  summarise(
    meanVal = mean(val, na.rm = T),
  ) %>%
  mutate(tn = as.numeric(tn)) %>%
  filter(ratio%in%c(0.1, 0.5, 0.9),
         stab%in%c(5, 10, 15),
         tn>1
         # measure%in%c("Mf")
  ) %>%
  ungroup()

ggplot(plotDat) + 
  facet_grid(~d, scales="free") +
  geom_path(aes(tn, meanVal, 
                colour=as.factor(ratio),
                linetype=as.factor(stab)),
            size=1) +
  my_theme
ggsave("../graphs/size_fitness2.pdf") 