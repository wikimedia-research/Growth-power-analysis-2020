## In this scenario, we're adding frwiki to the set of eight we currently have.

## Dataset of monthly averages of registrations, activations, and retention
## across 2019.
user_data = fread('datasets/aggregate_statistics.tsv')

## The group of wikis we're working with in this part of the analysis:
wikis = c('cswiki', 'kowiki', 'arwiki', 'viwiki',
          'euwiki', 'huwiki', 'ukwiki', 'hywiki',
          'frwiki')

## We'll round registration, activation, and retention down to the nearest integer.
user_data$n_registered = floor(user_data$n_registered)
user_data$n_activated = floor(user_data$n_activated)
user_data$n_retained = floor(user_data$n_retained)

## beta is calculated by reversing the "divide by 4 rule" from Gelman and Hill, so 10% x 4:
beta <- 0.1 * 4

## we choose a sigma_beta that's almost including 0 if multiplied by 1.96
sigma_beta <- 0.2

## The weighted average desktop activation rate from our wikis is the default.
## This gives the following mu:
mu <- logit(
  user_data[wiki_db %in% wikis & platform == 'desktop',
            list(avg_act = sum(n_activated)/sum(n_registered))]$avg_act
)

## We can estimate sigma_mu using the user data above, but in this case I want a wide range
## of possible values for mu, while also something that allows for variation in mobile.
## This seems like a reasonable compromise.
sigma_mu <- 0.3
arm::invlogit(mu + 3 * c(-sigma_mu, sigma_mu))

## We assume Normality for alphas, so we need a sigma for those:
sigma_alpha <- 0.2

# The 99.7% CI of alpha then becomes:
arm::invlogit(mu + 3 * c(-sigma_alpha, sigma_alpha))
## This CI doesn't cover the lowest value of desktop registration found in the dataset,
## but if we were to do so we'd also go way above the largest value. So we're making
## a sensible compromise.

## For the effect of mobile, we'll again calculate a weighted average:
user_data[wiki_db %in% wikis & platform == 'mobile',
          list(avg_act = sum(n_activated)/sum(n_registered))]$avg_act

## We then want (mu + mobile) to result in the given average probability:
beta_mob <- -0.11
arm::invlogit(mu + beta_mob)

# In the varying-intercept, varying-slope model the slope and intercept come from
# a MVN distribution and have between-group correlation. This is a messy & inexact
# way of guessing at that correlation for use in the model. See for more info:
# http://www.bristol.ac.uk/cmm/learning/videos/random-slopes.html#covar
prop_acts <- as.data.frame(
  user_data[wiki_db %in% wikis,
            list(prop_act = n_activated/n_registered),
            by = c("wiki_db", "platform")]
)

# First we plot the activation lines to see how is_mobile affects the raw proportion:
par(mfrow = c(1, 2))
plot(NA, NA, xlim = c(0, 1), ylim = c(0.1, 0.6), xlab = "platform", ylab = "prop_act",
     main = "Per-wiki differences of activation\nproportions by display")
wiki_colors <- brewer.pal(n = length(wikis), name = "Set1")
names(wiki_colors) = wikis  
for (wiki in wikis) {
  x0 <- 0
  y0 <- prop_acts$prop_act[prop_acts$wiki_db == wiki & prop_acts$platform == 'desktop']
  x1 <- 1
  y1 <- prop_acts$prop_act[prop_acts$wiki_db == wiki & prop_acts$platform == 'mobile']
  segments(x0, y0, x1, y1, col = wiki_colors[wiki], lwd = 2)
  points(c(x0, x1), c(y0, y1), col = wiki_colors[wiki], pch = 16, cex = 2)
  text(x0, y0 + 0.01, sprintf("%.1f%%", 100 * y0), col = wiki_colors[wiki], pos = 4, offset = 1)
  text(x1, y1 + 0.01, sprintf("%.1f%%", 100 * y1), col = wiki_colors[wiki], pos = 2, offset = 1)
  slope <- (y1 - y0)
  text(0.5, y0 + 0.5 * slope, sprintf("%.3f", slope), col = wiki_colors[wiki], pos = 3)
}; rm(x0, x1, y0, y1)

# Then we calculate a "slope" and an "intercept" and plot them to see if maybe there's
# a pattern. We don't have a lot of groups to work with so these results are iffy at
# best, but we can *try* to see if bigger intercept = bigger slope (positive correlation)
# or maybe bigger intercept = smaller slope (negative correlation):
slopes = c()
intercepts = c()

plot(NA, NA, xlim = c(-2, 0), ylim = c(-0.8, 0.8), xlab = "intercept", ylab = "slope",
     main = "Per-wiki intercept vs slope")
for (wiki in wikis) {
  x0 <- 0
  y0 <- prop_acts$prop_act[prop_acts$wiki_db == wiki & prop_acts$platform == 'desktop']
  x1 <- 1
  y1 <- prop_acts$prop_act[prop_acts$wiki_db == wiki & prop_acts$platform == 'mobile']
  intercept <- logit(y0)
  slope <- logit(y1) - logit(y0)
  slopes = c(slopes, slope)
  intercepts = c(intercepts, intercept)
  abline(h = 0, lty = "dashed")
  points(intercept, slope, col = wiki_colors[wiki], pch = 16, cex = 2)
}; rm(x0, x1, y0, y1, slope, intercept)
# Reset:
par(mfrow = c(1, 1))

## Possible correlation is similar to what we saw in the previous iteration ("current set"),
## so let's calculate a correlation both with and without the low slope value.
wiki_sl_int = data.table(slope = slopes, intercept = intercepts)
cor(wiki_sl_int$slope, wiki_sl_int$intercept)

## Correlation is reported as -0.6. I think that's slightly high due to
## the large value at the bottom, so let's remove that and try again.
cor(wiki_sl_int[slope > -0.5]$slope, wiki_sl_int[slope > -0.5]$intercept)

## It's at -0.6, but I want to be slightly conservative:
alpha_beta_rho <- -0.5

## We have both positive and negative effects of mobile, so sigma_mob should
## be chosen to reasonably simulate both.
sigma_mob <- 0.12
arm::invlogit(mu + (beta_mob + 3 * c(-sigma_mob, sigma_mob)))

## Lastly, we need a probability of registering on the mobile site, with variation on
## a per-wiki basis.
p_mob <- user_data[wiki_db %in% wikis,
                   list(prop_mobile = n_registered[platform == 'mobile'] / sum(n_registered)),
                   by = "wiki_db"]$prop_mobile



## Turn all our coefficients into a list:
act_coeff_addfr <- list(
  beta = beta,
  sigma_beta = sigma_beta,
  mu = mu,
  sigma_mu = sigma_mu,
  sigma_alpha = sigma_alpha,
  beta_mob = beta_mob,
  sigma_mob = sigma_mob,
  alpha_beta_rho = alpha_beta_rho,
  p_mob = p_mob
)

# Calculate once and reuse registration counts (note that we assume a 30-day month)
reg_counts_addfr <- user_data[wiki_db %in% wikis,
                              list(n_registered = sum(n_registered)),
                              by = wiki_db]$n_registered/30

save(list = c("act_coeff_addfr", "reg_counts_addfr"),
     file = "simulation_parameters/addfr_params.txt",
     ascii = TRUE)

# Wrapper for random numbers from bivariate normal:
rbvnorm <- function(n, mus, sigmas, corr) {
  mvtnorm::rmvnorm(
    n = n,
    mean = mus,
    sigma = matrix(
      c(sigmas[1]^2, rep(corr * sigmas[1] * sigmas[2], 2), sigmas[2]^2),
      nrow = 2, byrow = TRUE
    )
  )
}
# Test it with strong correlations:
plot(rbind(
  rbvnorm(500, c(mu, beta_mob), c(sigma_alpha, sigma_mob), 0.9),
  rbvnorm(500, c(mu, beta_mob), c(sigma_alpha, sigma_mob), -0.9)
), pch = 16, xlab = "intercept", ylab = "slope")

simulate_dataset <- function(N, sim_coeffs) {
  ## Simulate a dataset for the given number of wikis and participant counts (N)
  ## using the given set of coefficients (sim_coeffs)

  ## Overall intercept
  mu_d <- rnorm(1, sim_coeffs[['mu']], sim_coeffs[['sigma_mu']])

  alpha_beta <- rbvnorm(
    n = length(N),
    mus = c(mu_d, sim_coeffs[['beta_mob']]),
    sigmas = c(sim_coeffs[['sigma_alpha']], sim_coeffs[['sigma_mob']]),
    corr = sim_coeffs[['alpha_beta_rho']]
  )

  ## Change in overall intercept for each wiki:
  alphas <- alpha_beta[, 1]

  ## Effect of being on mobile for each wiki:
  beta_mobs <- alpha_beta[, 2]

  ## Treatment effect:
  beta_d <- rnorm(1, sim_coeffs[['beta']], sim_coeffs[['sigma_beta']])

  fake_data <- imap_dfr(N, ~ data.frame(
    wiki = .y,
    treatment = rbinom(.x, 1, 0.5), # 50/50 sampling into treatment vs control
    is_mobile = rbinom(.x, 1, p_mob[.y]) # per-wiki mobile registration probability
  )) %>% dplyr::mutate(
    p = invlogit(alphas[wiki] + beta_d * treatment + beta_mobs[wiki] * is_mobile),
    y = map_int(p, ~ rbinom(1, 1, prob = .x))
  )
  return(fake_data)
}

estimate_power <- function(N, sim_coeffs, n_sims = 1000) {
  print(paste("running", n_sims, "simulations for", length(N), "wikis, with", paste(N, collapse=', '), "users"))
  significant <- future_map_lgl(1:n_sims, function(i) {
    fit <- lme4::glmer(
      y ~ treatment + (1 + is_mobile | wiki),
      data = simulate_dataset(N, sim_coeffs),
      family = binomial()
    )
    beta_est <- fixef(fit)["treatment"]
    beta_se <- se.fixef(fit)["treatment"]
    return(beta_est - 2 * beta_se > 0) # H0: beta <= 0
  })
  return(mean(significant)) # Pr(reject H0 | H0 is false)
}


## We switch to calculating activation power on a weekly basis for up to 12 weeks,
## so our multiplier becomes the number of days in that.
multipliers <- seq(7, 84, by = 7)
# So the output is also named, with the estimated power for each multiplier:
multipliers <- setNames(multipliers, multipliers)
# Then use the ~ syntax to define an anonymous function applied to each multiplier:
activation_powers = list()
activation_powers[['10%']] <- purrr::map_dbl(
  multipliers,
  ~ estimate_power(floor(registration_counts * .x), act_coeff_first4)
)

## What if the effect is 5%?
activation_coefficients[['beta']] = 0.2
activation_coefficients[['sigma_beta']] = 0.1

activation_powers[['5%']] <- purrr::map_dbl(
  multipliers,
  ~ estimate_power(floor(registration_counts * .x), activation_coefficients)
)

## What if the effect is 2%?
activation_coefficients[['beta']] = 0.08
activation_coefficients[['sigma_beta']] = 0.04

activation_powers[['2%']] <- purrr::map_dbl(
  multipliers,
  ~ estimate_power(floor(registration_counts * .x), activation_coefficients)
)

## Flatten the list into a data.frame
x = rbindlist(lapply(activation_powers, as.data.frame.list))
x[, `Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))]

x %>% gather(Days, Power, -`Effect size`) %>%
  dplyr::mutate(num_days = as.numeric(stringr::str_extract(Days, "\\d+$"))) %>%
  ggplot(aes(
    x = factor(num_days), group = factor(`Effect size`),
    y = Power, color = factor(`Effect size`)
  )) +
  stat_identity(geom = "line", size = 1.1) +
  geom_point() +
  geom_label(aes(label = scales::percent(Power, 0.1)), show.legend = FALSE) +
  labs(
    x = "Number of days in experiment", y = "Power", color = "Effect size",
    title = "Target wikis Activation A/B Test power analysis"
  ) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(labels = scales::percent_format(1), limits = c(0, 1)) +
  theme(legend.position = "bottom")

