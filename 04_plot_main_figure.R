# PLOT MAIN FIGURE

# Load all model outputs and create multi-panel plot


# FINAL PLOTTING STYLE:
fit_no_na %>%
  spread_draws(tx_land_feed_w[tx]) %>%
  median_qi(tx_land_feed_w, .width = c(0.95, 0.8, 0.5)) %>%
  left_join(tx_index_key, by = "tx") %>% # Join with index key to get sci and taxa names
  # SET plants and bivalves to 0
  mutate(tx_land_feed_w = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = tx_land_feed_w),
         .lower = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .lower),
         .upper = if_else(taxa %in% c("bivalves", "plants"), true = 0, false = .upper)) %>%
  ggplot(aes(y = full_taxa_name, x = tx_land_feed_w)) +
  # can adjust intervals for thick vs. thin lines: stat_halfeye(.width = c(.80, .95))
  #stat_halfeye(aes(slab_fill = full_taxa_name), normalize = "all") +
  #stat_interval(aes(x = tx_land_feed_w), .width = c(0.5, 0.8, 0.95)) +
  geom_interval(aes(xmin = .lower, xmax = .upper)) +
  geom_point(x = c(0), y = c("bivalves")) +
  geom_point(x = 0, y = "plants") +
  scale_color_brewer() +
  theme_classic() + 
  tx_plot_theme + 
  labs(x = units_for_plot, y = "", title = "Mean off-farm (feed) impact")
ggsave(filename = file.path(outdir, paste("plot_", impact, "_", set_allocation, "-allocation_OFF-FARM-TAXA-LEVEL-WEIGHTED.png", sep = "")), width = 11, height = 8.5)
#ggsave(filename = file.path(outdir, "plot_geom_interval.png"))
