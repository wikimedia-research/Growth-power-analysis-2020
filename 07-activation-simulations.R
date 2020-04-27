## Run all simulations to calculate power for newcomer activation

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

  ## Probability of registering on mobile:
  p_mob = sim_coeffs[['p_mob']]
  
  fake_data <- purrr::imap_dfr(N, ~ data.frame(
    wiki = .y,
    treatment = rbinom(.x, 1, 0.5), # 50/50 sampling into treatment vs control
    is_mobile = rbinom(.x, 1, p_mob[.y]) # per-wiki mobile registration probability
  )) %>% dplyr::mutate(
    p = arm::invlogit(alphas[wiki] + beta_d * treatment + beta_mobs[wiki] * is_mobile),
    y = purrr::map_int(p, ~ rbinom(1, 1, prob = .x))
  )
  return(fake_data)
}

simulated_model = function(sim_n, N, sim_coeffs) {
  fit <- lme4::glmer(
      y ~ treatment + (1 + is_mobile | wiki),
      data = simulate_dataset(N, sim_coeffs),
      family = binomial()
  )
  beta_est <- lme4::fixef(fit)["treatment"]
  beta_se <- arm::se.fixef(fit)["treatment"]
  return(beta_est - 2 * beta_se > 0) # H0: beta <= 0
}

estimate_power <- function(N, sim_coeffs, n_sims = 1000) {
  print(paste("running", n_sims, "simulations for", length(N), "wikis, with",
              paste(N, collapse=', '), "users")
  )
  significant <- future_map_lgl(
    1:n_sims,
    simulated_model,
    N, sim_coeffs
  )
  return(mean(significant)) # Pr(reject H0 | H0 is false)
}

run_power_simulations = function(reg_counts, coeffs) {
  ## Run simulations of user activation with the given registraction counts
  ## and coefficients. Returns a flattened data frame with results for
  ## 10%, 5%, and 2% effect size.
  
  ## Activation power is calculated on weekly basis for up to 12 weeks.
  multipliers <- seq(7, 84, by = 7)
  # So the output is also named, with the estimated power for each multiplier:
  multipliers <- setNames(multipliers, multipliers)
  
  # Then use the ~ syntax to define an anonymous function applied to each multiplier:
  activation_powers = list()
  activation_powers[['10%']] <- purrr::map_dbl(
    multipliers,
    ~ estimate_power(floor(reg_counts * .x), coeffs)
  )
  
  ## What if the effect is 5%?
  coeffs[['beta']] = 0.2
  coeffs[['sigma_beta']] = 0.1
  
  activation_powers[['5%']] <- purrr::map_dbl(
    multipliers,
    ~ estimate_power(floor(reg_counts * .x), coeffs)
  )
  
  ## What if the effect is 2%?
  coeffs[['beta']] = 0.08
  coeffs[['sigma_beta']] = 0.04
  
  activation_powers[['2%']] <- purrr::map_dbl(
    multipliers,
    ~ estimate_power(floor(reg_counts * .x), coeffs)
  )
  
  ## Flatten the list into a data.frame
  x = rbindlist(lapply(activation_powers, as.data.frame.list))
  x[, `Effect size` := ordered(c("10%", "5%", "2%"), c("10%", "5%", "2%"))]
  return(x)
}

## 1: Load in the parameters
load('simulation_parameters/first4_params.txt', verbose = TRUE)
## 2: Run simulations
# first4_results = run_power_simulations(reg_counts_first4, act_coeff_first4)

## 3: save output as a TSV
# write.table(first4_results, file = 'datasets/first4-results.tsv', sep = '\t', quote = FALSE,
#            row.names = FALSE)

## Do this for all four configurations (would've been easier if the variable names
## were the same in all four now, wouldn't it?)

