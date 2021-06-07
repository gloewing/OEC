##
library(tidyverse)

##
eta <- read_csv("Manuscript Figures/Eta Figure/etaData") %>%
  # read_csv("Figures/Eta/etaData") %>%
  select(-X1) %>%
  pivot_longer(-c(eta, type), "week", "value") %>%
  mutate(week = as.numeric(str_replace_all(week, "week_", "")))
  
##
obs <- eta %>%
  filter(type == "obs")

##
preds <- eta %>%
  filter(type == "preds")

##
stack <- eta %>%
  filter(type == "stacking")

##
cspec <- eta %>%
  filter(type == "country_specific_model")

##
ggplot() +
  labs(x = "Weeks",
       y = "Rate of mortality") +
  geom_vline(xintercept = 52,
             lty = 2, 
             color = "gray") +
  geom_point(aes(week, value),
             alpha = 0.30, 
             data  = filter(obs, week <= 52)) +
  geom_point(aes(week, value),
             alpha = 0.60, 
             data  = filter(obs, week > 52)) +
  scale_y_continuous(breaks = seq(18, 28, 2)) +
  geom_line(aes(week, value, color = eta, group = eta),
            data = preds) +
  scale_color_viridis_c(name = expression(eta),
                        limits    = c(0.0, 1)) +
  geom_line(aes(week, value),
            color = "#ca0020",
            size  = 1,
            data  = stack) +
  geom_line(aes(week, value),
            color = "#525252",
            size  = 1,
            data  = cspec) +
  geom_label(aes(80, 19, label = "Country specific model"),
             fontface = "bold",
             color    = "#ca0020") +
  geom_label(aes(80, 18, label = "Stacking"),
             fontface = "bold",
             color    = "#525252") +
  theme_bw() +
  guides(color = guide_colorbar(title.position = "top", 
                                title.hjust = .5, 
                                barwidth = unit(.5, "lines"), 
                                barheight = unit(10, "lines")))


##
ggsave(filename = "eta-plot.pdf",
       height = 3,
       width  = 5)