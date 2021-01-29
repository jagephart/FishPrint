# Summarize results from posteriors

# Set data directories
datadir <- "/Volumes/jgephart/BFA Environment 2/Data"
outdir <- "/Volumes/jgephart/BFA Environment 2/Outputs"

#_______________________________________________________________________________________________________________________#
# Load packages and source functions
#_______________________________________________________________________________________________________________________#
library(tidyverse)
library(ggplot2)

#_______________________________________________________________________________________________________________________#
# Load results files and merge
#_______________________________________________________________________________________________________________________#

df <- read.csv(file.path(outdir, "PRIORS/GHG/summary_Global warming potential_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
df <- df %>%
  mutate(source = "on-farm", stressor = "GHG") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")

tmp <- read.csv(file.path(outdir, "PRIORS/GHG/summary_Global warming potential_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>% 
  mutate(source = "off-farm", stressor = "GHG") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp)

tmp <- read.csv(file.path(outdir, "PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "N") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Nitrogen/summary_Marine eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "N") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "P") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Phosphorus/summary_Freshwater eutrophication_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "P") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Water/summary_Water consumption_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Water") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Water/summary_Water consumption_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Water") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Land/summary_Land Use_Mass-allocation_ON-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "on-farm", stressor = "Land") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "on_farm", "lower_95" = ".lower", "upper_95" = ".upper")
df <- df %>%
  bind_rows(tmp) 

tmp <- read.csv(file.path(outdir, "PRIORS/Land/summary_Land Use_Mass-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.csv"))
tmp <- tmp %>%
  mutate(source = "off-farm", stressor = "Land") %>%
  select(taxa, full_taxa_name, stressor, source, "median" = "off_farm", "lower_95" = ".lower", "upper_95" = ".upper")
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
base_size <- 10
base_family <- "sans"

png("stressor_by_source.png", width = 8.5, height = 5, units = "in", res = 300)
ggplot(df, aes(x = median, y = taxa, fill = source)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#3FC1C9", "#57D182")) +
  labs(y = "", x = "") +
  facet_wrap(~stressor, nrow = 1, scales = "free") +
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
dev.off()


# Capture fishery results
df_capture <- read.csv(file.path(outdir, "PRIORS/Wild/summary_WILD-GHG-TAXA-LEVEL-WEIGHTED.csv"))

df_capture <- df_capture %>%
  filter(.width == 0.95)

df_total %>% 
  filter(stressor == "GHG") %>%
  arrange(total)


#_______________________________________________________________________________________________________________________#
# Write out results tables
#_______________________________________________________________________________________________________________________#
df_capture <- df_capture %>% 
  mutate(stressor = "GHG", full_taxa_name = taxa) %>%
  mutate(source = "All", Production = "Capture") %>%
  select(Production, taxa, full_taxa_name, stressor, source,  "median" = "total_stressor", 
         "lower_95" = ".lower", "upper_95" = ".upper")
  
df_all <- df %>%
  mutate(Production = "Aquaculture") %>%
  bind_rows(df_capture)

write.csv(df_all, file.path(outdir, "SI_stressor_results.csv"), row.names = FALSE)
