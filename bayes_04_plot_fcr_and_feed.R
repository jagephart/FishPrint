# PLOT other intermediate-level calculations (these are universally shared among all models)

# FIRST, load model output
# FCR and Feed will be the same for C, N, P, land, and water so just pick one
# But run separately for the different allocation methods: Mass, economic, gross energy content

# FCR
# Mean FCR sci-level
fit_no_na %>%
  spread_draws(sci_mu_fcr[sci]) %>%
  median_qi(.width = 0.95) %>%
  left_join(sci_index_key, by = "sci") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  mutate(sci_mu_fcr = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = sci_mu_fcr),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name))) %>%
  ggplot(aes(y = clean_sci_name, x = sci_mu_fcr, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
  geom_pointinterval() +
  theme_classic() + 
  sci_plot_theme + 
  labs(x = "FCR", y = "", title = "Mean FCR", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_FCR-SCI-LEVEL.png", sep = "")), width = 11, height = 8.5)


# Mean FCR taxa-level
fit_no_na %>%
  spread_draws(tx_mu_fcr[tx]) %>%
  median_qi(.width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  mutate(tx_mu_fcr = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_mu_fcr),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  # REORDER taxa axis
  mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order)) %>%
  #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_mu_fcr)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = c(0), y = c("bivalves")) +
  geom_point(x = 0, y = "plants") +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = "FCR", y = "", title = "Mean FCR", color = "taxa group")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_FCR-TAXA-LEVEL.png", sep = "")), width = 11, height = 8.5)

# FEED PROPORTIONS:
# For feed proportions, remove scinames where all studies were fcr == 0 and feed proportions were modified to 0.25
sci_null_feed <- lca_model_dat %>% 
  group_by(clean_sci_name) %>% 
  mutate(drop_feed_model = if_else(fcr == 0 & feed_soy == 0.25 & feed_crops == 0.25 & feed_fmfo == 0.25 & feed_animal == 0.25, true = 1, false = 0)) %>%
  mutate(n_sci = n()) %>%
  mutate(n_drop = sum(drop_feed_model)) %>%
  ungroup() %>%
  filter(n_drop == n_sci) %>% # the number of studies that were modified to 0.25 is the same as the total number of studies for that sci name
  pull(clean_sci_name) %>%
  unique()

# Same with taxa names
taxa_null_feed <- lca_model_dat %>% 
  group_by(taxa) %>% 
  mutate(drop_feed_model = if_else(fcr == 0 & feed_soy == 0.25 & feed_crops == 0.25 & feed_fmfo == 0.25 & feed_animal == 0.25, true = 1, false = 0)) %>%
  mutate(n_sci = n()) %>%
  mutate(n_drop = sum(drop_feed_model)) %>%
  ungroup() %>%
  filter(n_drop == n_sci) %>%
  pull(taxa) %>%
  unique()

# Sci-level theta (feed proportions)
# Make separate plot for each feed component
feed_component <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_component)){
  plot_dat <- fit_no_na %>%
    spread_draws(sci_theta[sci, feed_index]) %>%
    median_qi(.width = 0.95) %>%
    left_join(sci_feed_key, by = c("sci" = "sci", "feed_index" = "feed_index")) %>% # Join with index key to get sci and taxa names
    filter(feed == feed_component[i]) %>%
    filter(clean_sci_name %in% sci_null_feed == FALSE) %>%
    ungroup() %>% # need to remove automatic sci-level grouping for fct_reorder to work
    mutate(clean_sci_name = fct_reorder(clean_sci_name, as.character(full_taxa_name)))
  p <- ggplot(plot_dat, aes(y = clean_sci_name, x = sci_theta, xmin = .lower, xmax = .upper, color = full_taxa_name)) +
    geom_pointinterval() +
    #stat_halfeye(aes(slab_fill = full_taxa_name)) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_classic() + 
    sci_plot_theme + 
    labs(x = "Feed Proportion", y = "", title = feed_component[i], color = "taxa group")
  print(p)
  ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_FEED-", feed_component[i], "-SCI-LEVEL.png", sep = "")), width = 11, height = 8.5)
}

# Taxa-level theta (feed proportions)
# Make separate plot for each feed component
feed_component <- c("soy", "crops", "fmfo", "animal")
for (i in 1:length(feed_component)){
  plot_dat <- fit_no_na %>%
    spread_draws(tx_theta[tx, feed_index]) %>%
    median_qi(.width = c(0.95, 0.8, 0.5)) %>%
    left_join(tx_feed_key, by = c("tx" = "tx", "feed_index" = "feed_index")) %>% # Join with index key to get sci and taxa names
    filter(feed == feed_component[i]) %>%
    filter(taxa %in% taxa_null_feed == FALSE) %>%
    # REORDER taxa axis
    mutate(full_taxa_name = fct_relevel(full_taxa_name, full_taxa_name_order))
    #mutate(full_taxa_name = fct_reorder(full_taxa_name, tx_feed_fp_w)) %>%
  p <- ggplot(plot_dat, aes(y = full_taxa_name, x = tx_theta)) +
    geom_interval(aes(xmin = .lower, xmax = .upper)) +
    #stat_halfeye(aes(slab_fill = full_taxa_name)) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_classic() + 
    tx_plot_theme + 
    theme(legend.position = "none") +
    labs(x = "Feed Proportion", y = "", title = feed_component[i])
  print(p)
  ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_FEED-", feed_component[i], "-TAXA-LEVEL.png", sep = "")), width = 11, height = 8.5)
}
