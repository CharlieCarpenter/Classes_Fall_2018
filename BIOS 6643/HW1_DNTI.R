#####################
#
# Homework 1 
#
# Problems not to be turned in
#
#
#####################

library(tidyverse)
library(magrittr)

chol <- read.csv("~/Documents/Classes_Fall_2018/BIOS 6643/HW/HW1/cholesterol.csv")

chol %<>% mutate(diff = after - before)

change_mod <- lm(diff~1, data = chol)
summary(change_mod)

baseline_mod <- lm(after ~ before, data = chol)
summary(baseline_mod)

plot(resid(baseline_mod), main = "BAC mod")

plot(resid(change_mod), main = "CS mod")

hybrid_mod <- lm(diff ~ before, data = chol)
summary(hybrid_mod)
