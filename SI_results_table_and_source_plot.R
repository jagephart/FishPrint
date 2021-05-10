# Bind results from posteriors for SI

# Set data directories
rm(list=ls())
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#_______________________________________________________________________________________________________________________#
# Load packages and source functions
#_______________________________________________________________________________________________________________________#
library(tidyverse)
library(ggplot2)
library(ggpubr)

#_______________________________________________________________________________________________________________________#
# Load main results files and merge (Priors + Mass allocation + Edible Weight)
#_______________________________________________________________________________________________________________________#

df <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/GHG/summary_Global warming potential_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df <- df %>%
  mutate(source = "on-farm", stressor = "GHG", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/GHG/summary_Global warming potential_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp)

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Water/summary_Water consumption_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Water/summary_Water consumption_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land", weight = "edible", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df <- df %>%
  bind_rows(tmp) 

df_total <- df %>% 
  group_by(taxa, full_taxa_name, stressor) %>% 
  summarise(total = sum(median))
#_______________________________________________________________________________________________________________________#
# Summarize patterns for paper text
#_______________________________________________________________________________________________________________________#

# Lowest taxa across stressors
df_total %>% 
  group_by(stressor) %>%
  slice_min(total, n = 4)

df_total %>% 
  group_by(stressor) %>%
  slice_max(total, n = 3)

# Stressor correlations 
tmp <- df_total %>% 
  pivot_wider(names_from = stressor, values_from = total)
cor(tmp[,3:7])

# Percent of on- versus off-farm 
source_percent <- df %>%
  select(taxa, full_taxa_name, stressor, source, median) %>%
  pivot_wider(names_from = source, values_from = median) %>%
  mutate(total = `on-farm` + `off-farm`) %>%
  mutate(percent_onfarm = 100*(`on-farm`/total),
         percent_offfarm = 100*(`off-farm`/total))

# GHG results
source_percent %>% 
  filter(stressor == "GHG") %>%
  arrange(percent_onfarm)

source_percent %>% 
  filter(stressor == "GHG") %>%
  arrange(`on-farm`)

source_percent %>% 
  filter(stressor == "GHG") %>%
  arrange(total)

# Land results
source_percent %>% 
  filter(stressor == "Land") %>%
  arrange(`on-farm`)

source_percent %>% 
  filter(stressor == "Land") %>%
  arrange(`off-farm`)

source_percent %>% 
  filter(stressor == "Land") %>%
  arrange(percent_offfarm)

source_percent %>% 
  filter(stressor == "Land") %>%
  arrange(desc(total))

# Water results
source_percent %>% 
  filter(stressor == "Water") %>%
  arrange(percent_offfarm)

source_percent %>% 
  filter(stressor == "Water") %>%
  arrange(`off-farm`)

source_percent %>% 
  filter(stressor == "Water") %>%
  arrange(total)

# N and P results
source_percent %>% 
  filter(stressor == "N") %>%
  arrange(percent_offfarm)

source_percent %>% 
  filter(stressor == "P") %>%
  arrange(percent_offfarm)

source_percent %>% 
  filter(stressor == "N") %>%
  arrange(desc(total))

source_percent %>% 
  filter(stressor == "P") %>%
  arrange(desc(total))

# Stacked bar of impact by source
#facet_labs <- c("GHG kg CO2-eq", "Land", "N", "P", "Water")
#names(facet_labs) <- c("GHG", "Land", "N", "P", "Water")

df$taxa <- factor(df$taxa, levels = rev(c("misc_diad", "misc_marine", "shrimp", "milkfish",
                                      "tilapia", "catfish", "oth_carp", "trout", "salmon",
                                      "hypoph_carp", "plants", "bivalves")))

base_size <- 12
base_family <- "sans"

# GHG plot
ghg_plot <- ggplot(df %>% filter(stressor == "GHG"), 
                   aes(x = median, y = taxa, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  scale_y_discrete(labels = rev(c("misc diad", "misc marine", "shrimp", "milkfish",
                              "tilapia", "catfish", "misc carp", "trout", "salmon",
                              "silver/bighead", "seaweeds", "bivalves"))) +
  labs(title = "GHG", x = expression("kg CO"[2]~"t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_text(size = ceiling(base_size*0.7), colour = "black"),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# N plot
N_plot <- ggplot(df %>% filter(stressor == "N"), 
                   aes(x = median, y = taxa, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Nitrogen", x = expression("kg N-eq t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# P plot
P_plot <- ggplot(df %>% filter(stressor == "P"), 
                 aes(x = median, y = taxa, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Phosphorus", x = expression("kg P-eq t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# Land plot
land_plot <- ggplot(df %>% filter(stressor == "Land"), 
                 aes(x = median, y = taxa, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Land", x = expression("m"^2~"t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))

# Water plot
water_plot <- ggplot(df %>% filter(stressor == "Water"), 
                    aes(x = median, y = taxa, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#70468C", "#57D182")) +
  labs(title = "Water", x = expression("m"^"3"~"t"^"-1"), y = "") +
  theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"), 
        axis.text = element_blank(),
        axis.title = element_text(size = ceiling(base_size*0.8)), 
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        strip.background = element_rect(linetype = 1, fill = "white"), 
        strip.text = element_text(), 
        strip.text.x = element_text(vjust = 0.5), 
        strip.text.y = element_text(angle = -90), 
        strip.placement.y = "outside",
        legend.position="bottom",
        plot.title = element_text(size = ceiling(base_size*1.1), face = "bold"), 
        plot.subtitle = element_text(size = ceiling(base_size*1.05)))


# png(file.path(outdir, "plot_stressor_by_source.png"), width = 8.5, height = 5, units = "in", res = 300)
# ggarrange(ghg_plot, N_plot, P_plot, land_plot, water_plot, nrow = 1,
#           common.legend = TRUE, legend = "bottom")
# dev.off()


# Capture fishery results
df_capture <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Wild/summary_WILD-GHG-TAXA-LEVEL-WEIGHTED.csv"))

df_capture <- df_capture %>%
  filter(.width == 0.95)

df_total %>% 
  filter(stressor == "GHG") %>%
  arrange(total)

#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Mass allocation + Live Weight
#_______________________________________________________________________________________________________________________#

df_mass_live <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/GHG/summary_Global warming potential_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df_mass_live <- df_mass_live %>%
  mutate(source = "on-farm", stressor = "GHG", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/GHG/summary_Global warming potential_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp)

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Water/summary_Water consumption_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Water/summary_Water consumption_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land", weight = "live", allocation = "mass") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_mass_live <- df_mass_live %>%
  bind_rows(tmp) 

#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Economic allocation + Edible Weight
#_______________________________________________________________________________________________________________________#

df_economic_edible <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/GHG/summary_Global warming potential_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df_economic_edible <- df_economic_edible %>%
  mutate(source = "on-farm", stressor = "GHG", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/GHG/summary_Global warming potential_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp)

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Water/summary_Water consumption_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Water/summary_Water consumption_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Edible-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land", weight = "edible", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_edible <- df_economic_edible %>%
  bind_rows(tmp) 


#_______________________________________________________________________________________________________________________#
# Load SI results files and merge 
# Priors + Economic allocation + Live Weight
#_______________________________________________________________________________________________________________________#

df_economic_live <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/GHG/summary_Global warming potential_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df_economic_live <- df_economic_live %>%
  mutate(source = "on-farm", stressor = "GHG", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/GHG/summary_Global warming potential_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp)

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Nitrogen/summary_Marine eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Phosphorus/summary_Freshwater eutrophication_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Water/summary_Water consumption_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Water/summary_Water consumption_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Economic-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "Nature-submitted-2021-05/Bayesian-Means-Live-Weight/PRIORS/Land-FCR-Priors-Only/summary_Land Use_Economic-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land", weight = "live", allocation = "economic") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
df_economic_live <- df_economic_live %>%
  bind_rows(tmp) 


#_______________________________________________________________________________________________________________________#
# Write out results tables
#_______________________________________________________________________________________________________________________#
df_capture <- df_capture %>% 
  mutate(stressor = "GHG", full_taxa_name = taxa) %>%
  mutate(source = "all", production = "capture", weight = "live", allocation = "mass") %>%
  select(production, taxa, full_taxa_name, stressor, source,  "median" = "total_stressor", 
         "lower_95" = ".lower", "upper_95" = ".upper", weight, allocation)
  
df_all <- df %>%
  bind_rows(df_mass_live, df_economic_edible, df_economic_live) %>%
  mutate(production = "aquaculture") %>%
  bind_rows(df_capture)

write.csv(df_all, file.path(outdir, "SI_stressor_results.csv"), row.names = FALSE)
