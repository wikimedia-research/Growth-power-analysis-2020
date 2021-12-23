## Combining our simulation results from the three scenarios
## into a master table, with visualizations of each scenario
## based on a given effect size.

target_ret_res = fread('datasets/first4-retention-results.tsv')
current_ret_res = fread('datasets/current-retention-results.tsv')
addfr_ret_res = fread('datasets/addfr-retention-results.tsv')
addall_ret_res = fread('datasets/addall-retention-results.tsv')


## Flatten all lists and turn it into a master data.table
retention_power_scenarios = rbindlist(
  list(
    target_ret_res %>% 
      dplyr::mutate(`Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))) %>%
      gather(Days, Power, -`Effect size`) %>%
      dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+(\\.\\d+)?$")),
                    scenario := 'Target wikis'),
    current_ret_res %>% 
      dplyr::mutate(`Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))) %>%
      gather(Days, Power, -`Effect size`) %>%
      dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+(\\.\\d+)?$")),
                    scenario := 'Current wikis'),
    addfr_ret_res %>% 
      dplyr::mutate(`Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))) %>%
      gather(Days, Power, -`Effect size`) %>%
      dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+(\\.\\d+)?$")),
                    scenario := 'Add French'),
    addall_ret_res %>% 
      dplyr::mutate(`Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))) %>%
      gather(Days, Power, -`Effect size`) %>%
      dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+(\\.\\d+)?$")),
                    scenario := 'All 16 wikis')
  )
)
retention_power_scenarios[
  , scenario := ordered(scenario, c('Target wikis', 'Current wikis',
                                    'Add French', 'All 16 wikis'))]

for(eff_size in c('10%', '5%', '2%')) {
  g = ggplot(retention_power_scenarios[`Effect size` == eff_size],
             aes(
               x = factor(num_days), group = scenario,
               y = Power, color = scenario
             )) +
    stat_identity(geom = "line", size = 1.1) +
    geom_point() +
    geom_label(aes(label = scales::percent(Power, 0.1)), show.legend = FALSE) +
    labs(
      x = "Number of months in experiment", y = "Power", color = "scenario",
      title = paste("Growth Experiments Editor Retention Power Analysis -", eff_size, "Effect Size")
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(breaks = c(0:5)*0.2,
                       labels = scales::percent_format(1), limits = c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(paste0('graphs/retention_power_analysis_', sub('%', '', eff_size), "perc.png"), g,
         dpi = 150, width = 16, height = 9)
}; rm(eff_size, g)
