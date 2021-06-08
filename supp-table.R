##
library(tidyverse)
library(lubridate)
library(xtable)

##
dat <- read_csv("data/world-mort.csv")

##
tab <- dat %>%
  filter(lubridate::year(date) < 2020) %>%
  group_by(country) %>%
  summarize(`First year of training data` = round(min(lubridate::year(date))),
            `Years of training data`      = round(max(lubridate::year(date)) - min(lubridate::year(date)) + 1)) %>%
  ungroup() %>%
  arrange(desc(`Years of training data`))

##
xtable(tab)