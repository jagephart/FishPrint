# Bayesian analysis of marine mammal risk
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(rstan)
library(ggrepel)
#library(shinystan)
library(brms)
library(tidybayes)
library(plotrix) # std.error() function

# Mac
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

# Load data
mm_risk_raw <- read.csv(file.path(datadir,"marine_mammal_risk.csv"))
fui <- read.csv(file.path(datadir, "marine_mammal_fui.csv"))

# Join with fui data and tidy data
# FREQUENTIST VERSION:
mm_risk <- mm_risk_raw %>%
  filter(str_detect(Reference, pattern = "Brown")| str_detect(Reference, pattern = "Micheli")) %>%
  pivot_longer(high:low, names_to = "risk", values_to = "n.species") %>%
  filter(!(is.na(n.species))) %>%
  select("mm_species" = species, "species" = fui, gear, Region, risk, n.species) %>%
  left_join(fui, by = c("species", "gear")) %>%
  filter(!is.na(ghg)) #%>%

mm_risk_aveghg <- mm_risk %>%
  group_by(mm_species, gear, risk, n.species) %>%
  summarise(ghg.ave = mean(ghg), ghg.se = std.error(ghg, na.rm = TRUE))  %>%
  ungroup()

# Begin BAYESIAN version
# Create groups for Bayesian analysis
mm_risk_bayes  <- mm_risk %>%
  group_by(mm_species, gear) %>% # Only group by mm_species + gear
  mutate(group_index = cur_group_id())

# Set DATA
ghg <- mm_risk_bayes$ghg

# INDICES
N <- nrow(mm_risk_bayes)
N_GRP <- length(unique(mm_risk_bayes$group_index))
obs_to_grp <- mm_risk_bayes$group_index

# Set data for stan:
# NO PRIORS
stan_data <- list(N = N,
                  N_GRP = N_GRP,
                  obs_to_grp = obs_to_grp,
                  ghg = ghg)
## FIX IT - lower bound for grp_mu should be 0 (Or, switch to gamma distribution model)
# Normal distribution model:
stan_grouped <- 'data {
  int<lower=0> N;  // number of observations
  vector[N] ghg; // data
  int N_GRP; // number of groups
  int obs_to_grp[N]; // group indicators
}
parameters {
  real<lower=0> mu;
  real<lower=0> grp_sigma;
  vector<lower=0>[N_GRP] grp_mu;
  real<lower=0> sigma;
}

model {
  // likelihood
  grp_mu ~ normal(mu, sigma);
  ghg ~ normal(grp_mu[obs_to_grp], grp_sigma);
}'

mod <- stan_model(model_code = stan_grouped)
# Note: For Windows, apparently OK to ignore this warning message:
# Warning message:
#   In system(paste(CXX, ARGS), ignore.stdout = TRUE, ignore.stderr = TRUE) :
#   'C:/rtools40/usr/mingw_/bin/g++' not found

# Fit model:
# Set seed while testing
fit_mod <- sampling(object = mod, 
                    data = stan_data, 
                    cores = 4, 
                    seed = "11729", 
                    iter = 2500) 
                    #control = list(adapt_delta = 0.99, max_treedepth = 15))
#fit_no_na <- sampling(object = no_na_mod, data = stan_data, cores = 4, iter = 5000, control = list(adapt_delta = 0.99))
summary(fit_mod)$summary

# Use tidybayes + ggdist for finer control of aes mapping (instead of bayesplots) 
get_variables(fit_mod)

######################################################################################################
# PLOT final outputs

# Calculate biodiversity risk index
mm_riskindex <- mm_risk_aveghg %>%
  mutate(risk.index = ifelse(risk == "low", 1*n.species,
                             ifelse(risk == "medium", 2*n.species,
                                    ifelse(risk == "high", 3*n.species, NA)))) %>%
  group_by(mm_species, gear) %>%
  summarise(risk.index = sum(risk.index)) %>%
  ungroup()

# Use tidybayes to get credible intervals and join with marine mammal group name and risk info
mm_risk_ghg <- fit_mod %>%
  spread_draws(grp_mu[grp]) %>%
  median_qi(grp_mu, .width = c(0.95, 0.8, 0.5)) %>%
  left_join(mm_risk_bayes, by = c("grp" = "group_index")) %>% # Join with index key to get group names back
  left_join(mm_riskindex, by = c("mm_species", "gear")) # Join with risk index

# SET PLOT THEMES
# Set colors to match Figure 1 (see bayes_04_compile_posteriors.R)
interval_palette <- c("#9EA8B7", "#6A7A90", "#364F6B") # for testing
midwater_trawl_palette <- c("#FFB4C4", "#FF86A4", "#FC5185") # PINKS
gillnets_palette <- c("#FFA647", "#B57736", "#704B25") # ORANGES
bottom_trawl_palette <- c("#57D182", "#42955E", "#2E5C3C") # GREENS
traps_palette <- c("#3FC1C9", "#348A8F", "#275659") # LIGHT BLUES

mm_plot_theme <- list(theme(title = element_text(size = 6),
                            axis.title.x = element_text(size = 6),
                            axis.text=element_text(size=6, color = "black"),
                            legend.text = element_text(size = 6, color = "black"),
                            legend.position = "bottom",
                            legend.box.margin=margin(-10,-10,-10,-10)))
units_for_plot = bquote('kg'~CO[2]*' t'^-1)


# ggplot(mm_risk_ghg) +
#   geom_interval(aes(x = grp_mu, y = risk.index, xmin = .lower, xmax = .upper)) +
#   theme_classic() + 
#   mm_plot_theme + 
#   labs(x = units_for_plot, y = "Risk index", title = "") +
#   scale_color_manual(values = interval_palette)
#ggsave(filename = file.path(outdir, "plot_WILD-GHG-TAXA-LEVEL-WEIGHTED.png"), width = 11, height = 8.5)

# Try setting colors by specifying in a column of the dataframe
mm_risk_ghg_colors <- mm_risk_ghg %>%
  mutate(interval_color = case_when(.width == 0.95 & gear == "Midwater Trawls" ~ midwater_trawl_palette[1],
                                    .width == 0.80 & gear == "Midwater Trawls"~ midwater_trawl_palette[2],
                                    .width == 0.50 & gear == "Midwater Trawls" ~ midwater_trawl_palette[3],
                                    .width == 0.95 & gear == "Gillnets and Entangling Nets" ~ gillnets_palette[1],
                                    .width == 0.80 & gear == "Gillnets and Entangling Nets" ~ gillnets_palette[2],
                                    .width == 0.50 & gear == "Gillnets and Entangling Nets" ~ gillnets_palette[3],
                                    .width == 0.95 & gear == "Bottom Trawls" ~ bottom_trawl_palette[1],
                                    .width == 0.80 & gear == "Bottom Trawls" ~ bottom_trawl_palette[2],
                                    .width == 0.50 & gear == "Bottom Trawls" ~ bottom_trawl_palette[3],
                                    .width == 0.95 & gear == "Traps and Lift Nets" ~ traps_palette[1],
                                    .width == 0.80 & gear == "Traps and Lift Nets" ~ traps_palette[2],
                                    .width == 0.50 & gear == "Traps and Lift Nets" ~ traps_palette[3],
                                    .width == 0.95 ~ interval_palette[1],
                                    .width == 0.80 ~ interval_palette[2],
                                    .width == 0.50 ~ interval_palette[3],
                                    TRUE ~ "no color")) %>%
  mutate(interval_color = factor(interval_color, levels = c(midwater_trawl_palette, gillnets_palette, bottom_trawl_palette, traps_palette))) %>%
  select(-c(n.species, fui, ghg, risk)) %>%
  unique()

# Version with all color hues in legend
# FIX IT - Format legend labels with new lines and capitalization
# ggplot(mm_risk_ghg_colors) +
#   geom_interval(aes(x = grp_mu, y = risk.index, xmin = .lower, xmax = .upper, color = interval_color), show.legend = TRUE) +
#   geom_point(aes(x = grp_mu, y = risk.index)) +
#   geom_text(aes(x = grp_mu, y = risk.index, label = mm_species), hjust = 0.3, vjust = -1) +
#   theme_classic() +
#   mm_plot_theme +
#   labs(x = units_for_plot, y = "Risk index", title = "", color = "Gear") +
#   scale_color_manual(values = levels(mm_risk_ghg_colors$interval_color),
#                      labels = c("", "Midwater trawls", "",
#                                 "", "Gillnets and entangling nets", "",
#                                 "", "Bottom trawls", "",
#                                 "", "Traps and lift nets", "")) +
#   guides(color = guide_legend(keywidth = 0.1, keyheight = 0.1, default.unit = "inch"))

# geom_text_repel(data = . %>% mutate(label = if_else(plot_name_2 %in% interesting_points, true = plot_name_2, false = "")), 
#                 aes(label=label), box.padding = unit(0.5, "lines"), size = 4) +

# Version with simplified legend
png(file.path(outdir, "plot_Figure-3_bayes.png"), width = 89, height = 60, units = "mm", res = 300) # as per Nature formatting guidelines: 89 mm for single column; 183 mm for double column
ggplot(mm_risk_ghg_colors) +
  geom_interval(aes(x = grp_mu, y = risk.index, xmin = .lower, xmax = .upper, color = interval_color), show.legend = TRUE) +
  geom_point(aes(x = grp_mu, y = risk.index)) +
  #geom_text(aes(x = grp_mu, y = risk.index, label = mm_species), hjust = 0.3, vjust = -1, size = 2) +
  geom_text_repel(data = . %>% mutate(label = if_else(.width == 0.95, true = mm_species, false = "")), 
                   aes(x = grp_mu, y = risk.index, label = label), hjust = 0.3, vjust = -1.5, size = 2, segment.color = "transparent") +
  theme_classic() +
  mm_plot_theme +
  labs(x = units_for_plot, y = "Risk index", title = "", color = "Gear") +
  scale_color_manual(values = levels(mm_risk_ghg_colors$interval_color),
                     breaks = c("#FF86A4", "#B57736", "#42955E", "#348A8F"),
                     labels = c("Midwater\nTrawls", "Gillnets and\nEntangling\nNets", "Bottom \nTrawls", "Traps and\nLift Nets")) 
dev.off()




# BACK TO JG's code:
# Plots with combined biodiversity risk index
mm_riskindex_aveghg <- mm_risk_aveghg %>%
  mutate(risk.index = ifelse(risk == "low", 1*n.species,
                             ifelse(risk == "medium", 2*n.species,
                                    ifelse(risk == "high", 3*n.species, NA)))) %>%
  group_by(mm_species, gear, ghg.ave, ghg.se) %>%
  summarise(risk.index = sum(risk.index))

# ggplot(mm_riskindex_aveghg, aes(x = ghg.ave, y = risk.index)) +
#   geom_point() +
#   geom_errorbar(aes(xmin=ghg.ave-(1.96*ghg.se), xmax=ghg.ave+(1.96*ghg.se)), width=.1) +
#   labs(y = "Risk index", x = "kg CO2-eq per tonne", size = "N. Species", colour = "Risk") +
#   theme_clean() +
#   facet_wrap(~gear, nrow = length(unique(mm_riskindex_aveghg$gear)))

base_size <- 6
base_family <- "sans"

cols <- c("#57D182", "#FFD947", "#FFA647", "#70468C")


#png("Fig3.png", width = 89, height = 89, units = "mm", res = 300)
ggplot(mm_riskindex_aveghg, aes(x = ghg.ave, y = risk.index, colour = factor(gear))) +
  geom_point() +
  geom_errorbar(aes(xmin=ghg.ave-(1.96*ghg.se), xmax=ghg.ave+(1.96*ghg.se)), width=.1) +
  scale_color_manual(values = cols, labels = function(x) str_wrap(x, width = 12)) +
  geom_text(aes(x = ghg.ave, y = risk.index, label = mm_species), hjust = 0.2, vjust = -0.7, 
            size = 2, colour = "black") +
  labs(y = "Risk index", x = "kg CO2-eq per tonne", size = "N. Species", colour = "Risk") +
  guides(colour=guide_legend(nrow=1, byrow=TRUE)) + 
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_text(size = ceiling(base_size), colour = "black"),
        axis.title = element_text(size = ceiling(base_size)), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        strip.background = element_rect(linetype = 0), strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), strip.text.y = element_text(angle = -90), 
        legend.text = element_text(size = ceiling(0.9*base_size), family = "sans"), 
        legend.title = element_blank(), 
        legend.key = element_rect(fill = "white", colour = NA), 
        legend.position="bottom",
        legend.margin=margin(c(1,1,1,1)),
        plot.title = element_text(size = ceiling(base_size*1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1)))
#dev.off()

