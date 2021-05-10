library(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
my_theme <- theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

dat <- read_csv('../results/firstRun') %>% distinct(d, epsilon, gain, loss, .keep_all=TRUE)

datSimp <- dat %>% select(-ID, -tn, -tq) %>% 
  filter(d==0.5, epsilon==10, gain ==0.1, loss<=0.1) %>%
  gather("key", "val", starts_with("t", ignore.case=FALSE)) %>%
  extract(key, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)") %>%
  filter(tn<n) %>% 
  select(d, epsilon, gain, loss, measure, tq, tn, val) %>%
  spread(measure, val) %>%
  select(-Tr) %>%
  mutate(tYreal = tY*tF,
         tXreal = tX*tF,
         popX = tX*tF*as.numeric(tn),
         popY = tX*tF*as.numeric(tn),
         R = tR)

# datSimp <- dat %>%
#   gather("key", "val", starts_with("t", ignore.case=FALSE)) %>%
#   extract(key, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)") %>%
#   spread(measure, val) %>%
#   mutate(tYreal = tY*tF,
#          tXreal = tX*tF,
#          popX = tX*tF*as.numeric(tn),
#          popY = tX*tF*as.numeric(tn),
#          R = tR) %>%
#   gather("measure", "val", tF, tR, tW, tX, tY, tYreal, tXreal, popX, popY) %>%
#   # filter(tn==2) %>%
#   mutate(gain_loss = paste(gain, loss, sep="_"),
#          envRatio = gain/loss,
#          avgRat = 1/((1/gain + 1/loss)/2),
#          stability = 1/avgRat) %>%
#   mutate(treatFull = as.factor(paste(treat, gain_loss, sep="_"))) %>%
#   select(d,
#          epsilon,
#          treat,
#          gain_loss,
#          treatFull,
#          tq,
#          tn,
#          measure,
#          val,
#          envRatio,
#          stability,
#          R,
#          gain,
#          loss) %>%
#   ungroup()
# 
# dd <- datSimp %>% 
#   filter(d>0) %>%
#   group_by(measure, epsilon, d, gain, loss, tn) %>%
#   mutate(norm = val/sum(val), val = val) %>%
#   pivot_wider(names_from = measure, 
#               values_from = c("val", "norm")) %>%
#   mutate(
#     realX = norm_tF * val_tX, 
#     realY = norm_tF * val_tY,
#     popX = norm_tF * val_popX,
#     popY = norm_tF * val_popY
#   ) %>%
#   ungroup()
# small <- select(dat, -ID) %>% 
#   filter(d==0.5, epsilon==10, gain ==0.1, loss<=0.1) %>% 
#   distinct(d, epsilon, gain, loss, .keep_all=TRUE) %>% 
#   gather("key", "val", starts_with("t", ignore.case=FALSE)) %>%
#   extract(key, c("measure", "tq", "tn"), "(\\w{2})(\\d)(\\d\\b)") %>%
#   spread(measure, val) 

