# -- Set up 
library(tidyverse)
library(lubridate)
library(data.table)
# https://www.mortality.org/Public/STMF/Outputs/
# https://www.mortality.org/Public/STMF_DOC/STMFNote.pdf

# -- Loading country codes and names
names <- fread("https://www.mortality.org/Public/STMF/Outputs/countrynamelist.txt", header = FALSE) %>%
  as_tibble() %>%
  setNames(c("country", "countrycode"))

# -- Loading mortality data
dat <- fread("https://www.mortality.org/Public/STMF/Outputs/stmf.csv") %>% 
  as_tibble() %>%
  filter(Sex == "b") %>%
  select(CountryCode, Year, Week, DTotal, RTotal) %>%
  setNames(c("countrycode", "year", "week", "outcome", "rate")) %>%
  left_join(names, by = "countrycode") %>%
  mutate(population = 52 * round(outcome / rate))

# -- Vector of dates
dates <- seq(make_date(1990, 01, 01), make_date(2021, 12, 31), by = "days")

# -- Vector of countries
countries <- sort(unique(dat$country))

# -- Adding date column to dataset
dat <- tibble(date = dates) %>%
  mutate(week = week(date),
         year = year(date)) %>%
  right_join(dat, by = c("week", "year")) %>%
  group_by(week, year, country) %>%
  mutate(date = last(date)) %>%
  ungroup() %>%
  unique() %>%
  select(date, country, outcome, population) %>%
  arrange(country, date)

# -- Estimating population size via linear interpolation
dat <- map_df(countries, function(x){
  
  #
  tmp_pop <- round(approx(x    = filter(dat, country == x, week(date) == 26)$date,
                          y    = filter(dat, country == x, week(date) == 26)$population,
                          xout = filter(dat, country == x)$date,
                          rule = 2)$y)
  
  #
  filter(dat, country == x) %>%
    mutate(population = tmp_pop)
})
write.csv(dat, "~/Desktop/world-mort.csv", row.names = FALSE)
# dat <- read_csv("~/Downloads/world-mort.csv")
# dat %>%
#   ggplot(aes(date, population)) +
#   geom_line() +
#   facet_wrap(~country, scales = "free_y")


