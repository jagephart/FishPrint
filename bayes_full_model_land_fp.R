# Use Bayes gamma regression to estimate harvest and yield
# Calculate land footprint as harvest / yield

# First run: process_data_for_analysis.R, then clear environment other than:
rm(list=ls()[!(ls() %in% c("lca_dat_clean_groups", "datadir", "outdir"))])

# Standardize yield columns as tonnes per Ha
# 1000 kg in 1 tonne
# 10,000 sq meters in 1 hectare

# Check that no entries report both units:
# lca_dat_clean %>%
#   filter(is.na(Yield_t_per_Ha)==FALSE & is.na(Yield_kg_per_m3)==FALSE)

# ON HOLD: volume to area not possible without assumptions
# lca_dat_clean %>%
#   mutate(Yield_t_per_Ha == if_else(is.na(Yield_kg_per_m3)==FALSE, true = ))

