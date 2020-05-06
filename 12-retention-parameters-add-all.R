## Editor retention power analysis:
## This analysis is similar to 02-activation-power-target-wikis.R, except we are
## working on editor retention.
## In this analysis, we're assuming that editor activations stay constant, in the
## same way that we've previously assumed that registrations stay constant.

## Dataset of monthly averages of registrations, activations, and retention
## across 2019.
user_data = fread('datasets/aggregate_statistics.tsv')

## The group of wikis we're working with in this part of the analysis:
wikis = unique(user_data$wiki_db)

## We'll round registration, activation, and retention down to the nearest integer.
user_data$n_registered = floor(user_data$n_registered)
user_data$n_activated = floor(user_data$n_activated)
user_data$n_retained = floor(user_data$n_retained)

## beta is calculated by reversing the "divide by 4 rule" from Gelman and Hill, so 10% x 4:
beta <- 0.1 * 4

## we choose a sigma_beta that's almost including 0 if multiplied by 1.96
sigma_beta <- 0.2

## We'll use the weighted average desktop retention rate from our four target wikis as the default.
## This gives the following mu:
mu <- logit(
  user_data[wiki_db %in% wikis & platform == 'desktop',
            list(avg_ret = sum(n_retained)/sum(n_activated))]$avg_ret
)

user_data[wiki_db %in% wikis & platform == 'desktop',
          list(avg_ret = sum(n_retained)/sum(n_activated))]$avg_ret

## We can estimate sigma_mu using the user data above, but in this case I want a wide range
## of possible values for mu, while also something that allows for variation in mobile.
## This seems like a reasonable compromise.
sigma_mu <- 0.2
arm::invlogit(mu + 3 * c(-sigma_mu, sigma_mu))

## We assume Normality for alphas, so we need a sigma for those:
sigma_alpha <- 0.15

# The 99.7% CI of alpha then becomes:
arm::invlogit(mu + 3 * c(-sigma_alpha, sigma_alpha))

## Checking the data to see how it relates to the CI:
user_data[wiki_db %in% wikis & platform == 'desktop',
          list(ret_rate = n_retained/n_activated), by = wiki_db]

## This CI covers the the majority of the 2019 data. It does not capture
## the two smaller wikis that have a high retention rate. We'll have to
## live with that compromise.

## For the effect of mobile, we'll again calculate a weighted average:
user_data[wiki_db %in% wikis & platform == 'mobile',
          list(avg_ret = sum(n_retained)/sum(n_activated))]$avg_ret

## We then want (mu + mobile) to result in the given average probability:
beta_mob <- -0.5
arm::invlogit(mu + beta_mob)

# In the varying-intercept, varying-slope model the slope and intercept come from
# a MVN distribution and have between-group correlation. This is a messy & inexact
# way of guessing at that correlation for use in the model. See for more info:
# http://www.bristol.ac.uk/cmm/learning/videos/random-slopes.html#covar
prop_rets <- as.data.frame(
  user_data[wiki_db %in% wikis,
            list(prop_ret = n_retained/n_activated),
            by = c("wiki_db", "platform")]
)

# First we plot the activation lines to see how is_mobile affects the raw proportion:
par(mfrow = c(1, 2))
plot(NA, NA, xlim = c(0, 1), ylim = c(0.05, 0.4), xlab = "is_mobile", ylab = "prop_ret",
     main = "Per-wiki differences of retention\nproportions by display")
wiki_colors <- alphabet2(n = length(wikis))
names(wiki_colors) = wikis  
for (wiki in wikis) {
  x0 <- 0
  y0 <- prop_rets$prop_ret[prop_rets$wiki_db == wiki & prop_rets$platform == 'desktop']
  x1 <- 1
  y1 <- prop_rets$prop_ret[prop_rets$wiki_db == wiki & prop_rets$platform == 'mobile']
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

plot(NA, NA, xlim = c(-2.5, 0), ylim = c(-1.5, 0.75), xlab = "intercept", ylab = "slope",
     main = "Per-wiki intercept vs slope")
for (wiki in wikis) {
  x0 <- 0
  y0 <- prop_rets$prop_ret[prop_rets$wiki_db == wiki & prop_rets$platform == 'desktop']
  x1 <- 1
  y1 <- prop_rets$prop_ret[prop_rets$wiki_db == wiki & prop_rets$platform == 'mobile']
  intercept <- logit(y0)
  slope <- logit(y1) - logit(y0)
  slopes = c(slopes, slope)
  intercepts = c(intercepts, intercept)
  abline(h = 0, lty = "dashed")
  points(intercept, slope, col = wiki_colors[wiki], pch = 16, cex = 2)
}; rm(x0, x1, y0, y1, slope, intercept, wiki)
# Reset:
par(mfrow = c(1, 1))

## In this case, it seems we have a fairly strong correlation for most of
## the wikis. I still want to keep the two smaller wikis with exceptional
## retention out as outliers.

wiki_sl_int = data.table(slope = slopes, intercept = intercepts)

cor(wiki_sl_int$slope, wiki_sl_int$intercept)

cor(wiki_sl_int[intercept < -1.0]$slope, wiki_sl_int[intercept < -1.0]$intercept)

## We find a medium correlation of -0.27. I'll choose to be slightly
## conservative here, though, and set it to -0.2
alpha_beta_rho <- -0.2

## We see that most of the wikis have negative association between mobile
## and retention, but some are positive. This makes us choose a sigma_mob that can
## model both a positive and negative impact.
## The chosen value has a slightly higher upside than the data, and a less severe
## downside than seen in the data. But, we have to compromise, so we do.
sigma_mob <- 0.15
arm::invlogit(mu + (beta_mob + 3 * c(-sigma_mob, sigma_mob)))

## Lastly, we need a probability of registering on the mobile site, with variation on
## a per-wiki basis.
p_mob <- user_data[wiki_db %in% wikis,
                   list(prop_mobile = n_registered[platform == 'mobile'] / sum(n_registered)),
                   by = "wiki_db"]$prop_mobile


## Turn all our coefficients into a list:
ret_coeff_addall <- list(
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

# Calculate once and reuse activation counts (we work on a monthly basis here):
act_counts_addall <- user_data[wiki_db %in% wikis,
                               list(n_activated = sum(n_activated)),
                               by = wiki_db]$n_activated

save(list = c("ret_coeff_addall", "act_counts_addall"),
     file = "simulation_parameters/addall_retention_params.txt",
     ascii = TRUE)
