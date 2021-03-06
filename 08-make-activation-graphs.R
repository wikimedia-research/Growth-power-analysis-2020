## Load in the four datasets, combine them, make graphs.

target_act_res = fread('datasets/first4-results.tsv')
current_act_res = fread('datasets/current-results.tsv')
addfr_act_res = fread('datasets/addfr-results.tsv')
addall_act_res = fread('datasets/addall-results.tsv')

## Flatten all lists and turn it into a master data.table
activation_power_scenarios = rbindlist(
  list(
    target_act_res %>% 
      dplyr::mutate(`Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))) %>%
      gather(Days, Power, -`Effect size`) %>%
      dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+$")),
                    scenario := 'Target wikis'),
    current_act_res %>% 
      dplyr::mutate(`Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))) %>%
      gather(Days, Power, -`Effect size`) %>%
      dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+$")),
                    scenario := 'Current wikis'),
    addfr_act_res %>% 
      dplyr::mutate(`Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))) %>%
      gather(Days, Power, -`Effect size`) %>%
      dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+$")),
                    scenario := 'Add French'),
    addall_act_res %>% 
      dplyr::mutate(`Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))) %>%
      gather(Days, Power, -`Effect size`) %>%
      dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+$")),
                    scenario := 'All 16 wikis')
  )
)
activation_power_scenarios[
  , scenario := ordered(scenario, c('Target wikis', 'Current wikis',
                                    'Add French', 'All 16 wikis'))]

for(eff_size in c('10%', '5%', '2%')) {
  g = ggplot(activation_power_scenarios[`Effect size` == eff_size],
             aes(
               x = factor(num_days), group = scenario,
               y = Power, color = scenario
             )) +
    stat_identity(geom = "line", size = 1.1) +
    geom_point() +
    geom_label(aes(label = scales::percent(Power, 0.1)), show.legend = FALSE) +
    labs(
      x = "Number of days in experiment", y = "Power", color = "scenario",
      title = paste("Growth Experiments Editor Activation Power Analysis -", eff_size, "Effect Size")
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_y_continuous(breaks = c(0:5)*0.2,
                       labels = scales::percent_format(1), limits = c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(paste0('graphs/activation_power_analysis_', sub('%', '', eff_size), "perc.png"), g,
         dpi = 150, width = 16, height = 9)
}; rm(eff_size, g)
