library(rstan)
library(here)
here::i_am("figure_stanfit/check_competition_classifications.R")

load(here::here("stanfit", "stanfit_default_all_refactor.rda"))
load(here::here("stanfit", "stanfit_default_all_narrow.rda"))
load(here::here("stanfit", "stanfit_default_all_wide.rda"))

fits <- list(
  default = stanfit_default_all_refactor,
  narrow = stanfit_default_all_narrow,
  wide = stanfit_default_all_wide
)
species_pairs <- list(c(1, 2), c(1, 3), c(2, 3))

for (prior in names(fits)) {
  draws <- rstan::extract(
    fits[[prior]], pars = c("lambda", "alpha", "rho", "fdiff")
  )

  growth <- draws$lambda - 1
  rho <- draws$rho
  inverse_rho <- 1 / rho
  fitness <- pmax(draws$fdiff, 1 / draws$fdiff)

  terms <- list(
    "lambda - 1" = growth, alpha = draws$alpha, rho = rho,
    fdiff = draws$fdiff, "1 / rho" = inverse_rho, "folded fdiff" = fitness
  )
  
  for (term in names(terms)) {
    values <- terms[[term]]
    if (!all(is.finite(values) & values > 0)) {
      stop(sprintf("%s: non-finite or non-positive %s", prior, term), call. = FALSE)
    }
  }

  original <- ifelse(
    rho < fitness & fitness < inverse_rho, "coexistence",
    ifelse(rho > fitness & fitness > inverse_rho, "priority effect", "competitive exclusion")
  )

  log_growth <- log(growth)
  log_alpha <- log(draws$alpha)
  log_rho <- log_fdiff <- array(NA_real_, dim = dim(rho))
  for (pair in seq_along(species_pairs)) {
    i <- species_pairs[[pair]][1]
    j <- species_pairs[[pair]][2]
    alpha_ii <- log_alpha[, i, i, ]
    alpha_ij <- log_alpha[, i, j, ]
    alpha_ji <- log_alpha[, j, i, ]
    alpha_jj <- log_alpha[, j, j, ]

    log_rho[, pair, ] <- (alpha_ij + alpha_ji - alpha_ii - alpha_jj) / 2
    log_fdiff[, pair, ] <- log_growth[, j, ] - log_growth[, i, ] +
      (alpha_ij + alpha_ii - alpha_ji - alpha_jj) / 2
  }
  
  if (!all(is.finite(log_rho) & is.finite(log_fdiff))) {
    stop(sprintf("%s: non-finite log-scale quantities", prior))
  }

  log_fitness <- abs(log_fdiff)
  recomputed <- ifelse(
    log_fitness < -log_rho, "coexistence",
    ifelse(log_fitness < log_rho, "priority effect", "competitive exclusion")
  )

  mismatches <- which(original != recomputed, arr.ind = TRUE)
  if (nrow(mismatches) > 0) {
    stop(sprintf(
      "%s: %d mismatches. First mismatch at draw %d, pair %d, group %d.",
      prior, nrow(mismatches), mismatches[1, 1], mismatches[1, 2], mismatches[1, 3]
    ))
  }

  message(sprintf(
    "%s: all %s classifications matched (%d draws x %d pairs x %d groups).",
    prior, format(length(original), big.mark = ","),
    dim(original)[1], dim(original)[2], dim(original)[3]
  ))
}
