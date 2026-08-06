# Scenario and outcome-model definitions used in the manuscript.
# Source R/siga_core.R before sourcing this file.

calibrate_logistic_intercept <- function(patterns,
                                         profile_prob,
                                         coefficients,
                                         interaction = 0,
                                         target_risk) {
  patterns <- check_binary_factor_matrix(patterns)
  profile_prob <- as.numeric(profile_prob)
  coefficients <- as.numeric(coefficients)
  if (length(coefficients) != ncol(patterns)) {
    stop("coefficient length does not match factor count.", call. = FALSE)
  }
  linear_without_intercept <- drop(patterns %*% coefficients)
  if (ncol(patterns) >= 2L) {
    linear_without_intercept <- linear_without_intercept +
      interaction * patterns[, 1L] * patterns[, 2L]
  }
  objective <- function(alpha) {
    sum(profile_prob * plogis(alpha + linear_without_intercept)) - target_risk
  }
  uniroot(objective, interval = c(-30, 30), tol = 1e-13)$root
}

calibrate_common_log_odds <- function(patterns,
                                      profile_prob,
                                      alpha0,
                                      coefficients,
                                      interaction,
                                      target_risk_difference) {
  patterns <- check_binary_factor_matrix(patterns)
  coefficients <- as.numeric(coefficients)
  eta0 <- alpha0 + drop(patterns %*% coefficients)
  if (ncol(patterns) >= 2L) {
    eta0 <- eta0 + interaction * patterns[, 1L] * patterns[, 2L]
  }
  p0 <- plogis(eta0)
  objective <- function(theta) {
    sum(profile_prob * (plogis(eta0 + theta) - p0)) - target_risk_difference
  }
  theta <- uniroot(objective, interval = c(-30, 30), tol = 1e-13)$root
  list(alpha0 = alpha0, theta = theta, p0 = p0, p1 = plogis(eta0 + theta))
}

make_independent_design <- function(factor_count,
                                    n_per_group,
                                    pbc = 0.80,
                                    factor_prob_master = c(0.50, 0.40, 0.30, 0.20, 0.10)) {
  K <- as.integer(factor_count)
  patterns <- all_binary_patterns(K)
  factor_prob <- factor_prob_master[seq_len(K)]
  profile_prob <- profile_prob_independent(factor_prob, patterns)
  list(
    K = K,
    J = nrow(patterns),
    total_n = 2L * as.integer(n_per_group),
    n_per_group = as.integer(n_per_group),
    pbc = pbc,
    weights = rep(1, K + 1L),
    patterns = patterns,
    factor_prob = factor_prob,
    profile_prob = profile_prob,
    profile_type = "independent"
  )
}

make_full_grid_designs <- function() {
  designs <- list(
    make_independent_design(2L, 100L),
    make_independent_design(2L, 500L),
    make_independent_design(5L, 200L),
    make_independent_design(5L, 1000L)
  )
  names(designs) <- c("K2_n100", "K2_n500", "K5_n200", "K5_n1000")
  for (i in seq_along(designs)) {
    designs[[i]]$design_id <- i
    designs[[i]]$design_label <- names(designs)[i]
  }
  designs
}

continuous_power_effects <- data.frame(
  factor_count = c(2L, 2L, 5L, 5L),
  n_per_group = c(100L, 500L, 200L, 1000L),
  superiority = c(0.430, 0.190, 0.300, 0.135),
  noninferiority = c(0.230, -0.010, 0.100, -0.065),
  equivalence = c(0.000, 0.280, 0.180, 0.330),
  stringsAsFactors = FALSE
)

binary_power_effects <- data.frame(
  factor_count = c(2L, 2L, 5L, 5L),
  n_per_group = c(100L, 500L, 200L, 1000L),
  superiority = c(0.200, 0.100, 0.140, 0.070),
  noninferiority = c(0.100, 0.000, 0.050, -0.030),
  equivalence_lower = c(-0.210, -0.100, -0.150, -0.100),
  equivalence_upper = c(0.210, 0.100, 0.150, 0.100),
  equivalence = c(0.000, 0.000, 0.000, 0.040),
  stringsAsFactors = FALSE
)

# Seven scenario roles are generated for each design and outcome type.
make_full_grid_scenarios <- function() {
  designs <- make_full_grid_designs()
  rows <- list()
  scenario_id <- 0L
  for (outcome in c("continuous", "binary")) {
    for (design_name in names(designs)) {
      design <- designs[[design_name]]
      if (outcome == "continuous") {
        effects <- continuous_power_effects[
          continuous_power_effects$factor_count == design$K &
            continuous_power_effects$n_per_group == design$n_per_group, , drop = FALSE
        ]
        L <- -0.20
        U <- 0.45
        lower_eq <- -0.45
        upper_eq <- 0.45
      } else {
        effects <- binary_power_effects[
          binary_power_effects$factor_count == design$K &
            binary_power_effects$n_per_group == design$n_per_group, , drop = FALSE
        ]
        L <- -0.10
        lower_eq <- effects$equivalence_lower
        upper_eq <- effects$equivalence_upper
        U <- upper_eq
      }
      roles <- list(
        list(objective = "Superiority", measure = "Type I error",
             effect = 0, lower = 0, upper = NA_real_, alpha = 0.05),
        list(objective = "Superiority", measure = "Power",
             effect = effects$superiority, lower = 0, upper = NA_real_, alpha = 0.05),
        list(objective = "Non-inferiority", measure = "Type I error",
             effect = L, lower = L, upper = NA_real_, alpha = 0.025),
        list(objective = "Non-inferiority", measure = "Power",
             effect = effects$noninferiority, lower = L, upper = NA_real_, alpha = 0.025),
        list(objective = "Equivalence", measure = "Error at lower limit",
             effect = lower_eq, lower = lower_eq, upper = upper_eq, alpha = 0.05),
        list(objective = "Equivalence", measure = "Error at upper limit",
             effect = upper_eq, lower = lower_eq, upper = upper_eq, alpha = 0.05),
        list(objective = "Equivalence", measure = "Power",
             effect = effects$equivalence, lower = lower_eq, upper = upper_eq, alpha = 0.05)
      )
      for (role in roles) {
        scenario_id <- scenario_id + 1L
        rows[[scenario_id]] <- data.frame(
          scenario_id = scenario_id,
          outcome = outcome,
          design_name = design_name,
          factor_count = design$K,
          n_per_group = design$n_per_group,
          total_n = design$total_n,
          objective = role$objective,
          measure = role$measure,
          effect = as.numeric(role$effect),
          lower_boundary = as.numeric(role$lower),
          upper_boundary = as.numeric(role$upper),
          alpha = as.numeric(role$alpha),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

make_continuous_benchmark_model <- function(design, effect) {
  K <- design$K
  beta <- design$factor_prob
  interaction <- 0.25
  centered_mean <- drop(sweep(design$patterns, 2L, design$factor_prob, "-") %*% beta)
  centered_mean <- centered_mean + interaction *
    (design$patterns[, 1L] * design$patterns[, 2L] -
       design$factor_prob[1L] * design$factor_prob[2L])
  list(
    type = "continuous_benchmark",
    effect = as.numeric(effect),
    mu0 = centered_mean,
    outcome_sd = 1,
    d_at = function(boundary) rep(as.numeric(effect) - boundary, design$J)
  )
}

make_binary_benchmark_model <- function(design, effect) {
  K <- design$K
  coefficients <- c(0.35, -0.25, 0.20, -0.15, 0.10)[seq_len(K)]
  interaction <- 0.20
  alpha0 <- calibrate_logistic_intercept(
    design$patterns, design$profile_prob, coefficients, interaction, target_risk = 0.60
  )
  fitted <- calibrate_common_log_odds(
    design$patterns, design$profile_prob, alpha0, coefficients, interaction,
    target_risk_difference = as.numeric(effect)
  )
  list(
    type = "binary_benchmark",
    effect = as.numeric(effect),
    alpha0 = alpha0,
    theta = fitted$theta,
    p0 = fitted$p0,
    p1 = fitted$p1,
    d_at = function(boundary) fitted$p1 - fitted$p0 - boundary
  )
}

generate_profile_sequence <- function(n, patterns, profile_prob) {
  id <- sample.int(nrow(patterns), size = n, replace = TRUE, prob = profile_prob)
  list(id = id, X = patterns[id, , drop = FALSE])
}

generate_benchmark_trial <- function(design, model, seed_profile, seed_allocation, seed_outcome) {
  set.seed(normalise_seed(seed_profile))
  profile <- generate_profile_sequence(design$total_n, design$patterns, design$profile_prob)
  z <- absolute_minimization_assign(
    profile$X, pbc = design$pbc, weights = design$weights, seed = seed_allocation
  )
  A <- as.integer((z + 1L) / 2L)
  set.seed(normalise_seed(seed_outcome))
  if (model$type == "continuous_benchmark") {
    y0 <- model$mu0[profile$id] + rnorm(design$total_n, sd = model$outcome_sd)
    y <- y0 + A * model$effect
  } else if (model$type == "binary_benchmark") {
    probability <- ifelse(A == 1L, model$p1[profile$id], model$p0[profile$id])
    y <- rbinom(design$total_n, 1L, probability)
  } else {
    stop("Unsupported benchmark model type.", call. = FALSE)
  }
  list(X = profile$X, profile_id = profile$id, z = z, A = A, y = y)
}

scale_direction <- function(direction, target_max_abs = 1) {
  direction <- as.numeric(direction)
  maximum <- max(abs(direction))
  if (!is.finite(maximum) || maximum <= 0) return(rep(0, length(direction)))
  direction * target_max_abs / maximum
}

weighted_center <- function(x, weight) {
  x <- as.numeric(x)
  weight <- as.numeric(weight)
  x - sum(weight * x) / sum(weight)
}

scale_direction_for_binary_delta <- function(direction,
                                             boundary,
                                             target_max_abs_d = 0.80,
                                             max_abs_delta = 0.90,
                                             safety = 0.98) {
  direction <- scale_direction(direction, 1)
  candidate <- Inf
  positive <- direction > 0
  negative <- direction < 0
  if (any(positive)) {
    candidate <- min(candidate, min((max_abs_delta - boundary) / direction[positive]))
  }
  if (any(negative)) {
    candidate <- min(candidate, min((-max_abs_delta - boundary) / direction[negative]))
  }
  if (!is.finite(candidate) || candidate <= 0) {
    stop("No positive binary heterogeneity scale satisfies the probability bounds.", call. = FALSE)
  }
  scale <- min(target_max_abs_d, safety * candidate)
  d <- scale * direction
  delta <- boundary + d
  if (any(abs(delta) >= 1)) stop("Binary risk differences are outside (-1,1).", call. = FALSE)
  d
}

make_pair_path_supplemental_designs <- function() {
  master <- c(0.50, 0.40, 0.30, 0.20, 0.10)
  specification <- data.frame(
    design_id = 1:6,
    label = paste0("D", 1:6),
    K = c(2L, 2L, 5L, 5L, 2L, 2L),
    target_per_group = c(100L, 500L, 200L, 1000L, 100L, 500L),
    total_n = c(200L, 1000L, 400L, 2000L, 200L, 1000L),
    pbc = c(0.80, 0.80, 0.80, 0.80, 0.95, 0.95),
    profile_type = "independent",
    latent_strength = 0,
    stringsAsFactors = FALSE
  )
  out <- vector("list", nrow(specification))
  for (i in seq_len(nrow(specification))) {
    row <- specification[i, ]
    patterns <- all_binary_patterns(row$K)
    factor_prob <- master[seq_len(row$K)]
    profile_prob <- profile_prob_independent(factor_prob, patterns)
    out[[i]] <- list(
      design_id = row$design_id, label = row$label, K = row$K,
      J = nrow(patterns), target_per_group = row$target_per_group,
      total_n = row$total_n, pbc = row$pbc,
      weights = rep(1, row$K + 1L), profile_type = row$profile_type,
      latent_strength = row$latent_strength, patterns = patterns,
      factor_prob = factor_prob, profile_prob = profile_prob,
      correlation = profile_correlations(patterns, profile_prob)
    )
  }
  names(out) <- specification$label
  out
}

make_pair_path_supplemental_scenarios <- function() {
  out <- data.frame(
    scenario_id = 1:12,
    scenario_code = sprintf("S%02d", 1:12),
    set = "supplemental",
    design = rep(paste0("D", 1:6), each = 2L),
    outcome = rep(c("continuous", "binary"), 6L),
    model = c(
      "realistic", "common_log_odds", "realistic", "common_log_odds",
      "realistic", "common_log_odds", "realistic", "common_log_odds",
      "pure_symmetric", "symmetric_stress", "pure_symmetric", "symmetric_stress"
    ),
    direction = c(
      rep(c("first_factor_contrast", "common_log_odds"), 4L),
      rep("robust_max_ratio", 4L)
    ),
    scenario_class = c(rep("practical_heterogeneity", 8L), rep("strong_pair_path", 4L)),
    stringsAsFactors = FALSE
  )
  out$boundary <- ifelse(out$outcome == "continuous", -0.20, -0.10)
  out$alpha <- 0.025
  out
}

make_pair_path_supplemental_model <- function(scenario, design, calibration, direction_calibration = NULL) {
  first_factor <- weighted_center(design$patterns[, 1L], design$profile_prob)
  first_factor <- scale_direction(first_factor, 1)
  direction <- switch(
    scenario$direction,
    first_factor_contrast = first_factor,
    robust_max_ratio = {
      if (is.null(direction_calibration)) stop("An independent direction calibration is required.", call. = FALSE)
      generalized_pair_directions(direction_calibration$psi, design$profile_prob)$max_direction
    },
    common_log_odds = rep(0, design$J),
    stop("Unknown supplemental direction.", call. = FALSE)
  )
  boundary <- as.numeric(scenario$boundary)
  if (scenario$outcome == "continuous" && scenario$model == "pure_symmetric") {
    d <- weighted_center(scale_direction(direction, 1.00), design$profile_prob)
    return(list(type = "continuous_symmetric", boundary = boundary,
                delta = boundary + d, d = d, shared_noise_sd = 0.05))
  }
  if (scenario$outcome == "continuous" && scenario$model == "realistic") {
    d <- weighted_center(scale_direction(direction, 0.50), design$profile_prob)
    beta <- design$factor_prob
    mu0 <- drop(sweep(design$patterns, 2L, design$factor_prob, "-") %*% beta)
    mu0 <- mu0 + 0.25 *
      (design$patterns[, 1L] * design$patterns[, 2L] -
         design$factor_prob[1L] * design$factor_prob[2L])
    return(list(type = "continuous_realistic", boundary = boundary,
                delta = boundary + d, d = d, mu0 = mu0,
                outcome_sd = 1, individual_effect_sd = 0.25))
  }
  if (scenario$outcome == "binary" && scenario$model == "symmetric_stress") {
    d <- scale_direction_for_binary_delta(direction, boundary, 0.80, 0.90, 0.98)
    d <- weighted_center(d, design$profile_prob)
    delta <- boundary + d
    p0 <- (1 - delta) / 2
    p1 <- (1 + delta) / 2
    return(list(type = "binary_direct_probability", boundary = boundary,
                p0 = p0, p1 = p1, d = p1 - p0 - boundary))
  }
  if (scenario$outcome == "binary" && scenario$model == "common_log_odds") {
    coefficients <- c(0.35, -0.25, 0.20, -0.15, 0.10)[seq_len(design$K)]
    alpha0 <- calibrate_logistic_intercept(
      design$patterns, design$profile_prob, coefficients, 0.20, 0.60
    )
    fitted <- calibrate_common_log_odds(
      design$patterns, design$profile_prob, alpha0, coefficients, 0.20,
      target_risk_difference = boundary
    )
    return(list(type = "binary_common_log_odds", boundary = boundary,
                p0 = fitted$p0, p1 = fitted$p1,
                d = fitted$p1 - fitted$p0 - boundary))
  }
  stop("Unknown supplemental model.", call. = FALSE)
}

generate_pair_path_supplemental_trial <- function(design, model, seed_profile, seed_allocation, seed_outcome) {
  set.seed(normalise_seed(seed_profile))
  profile <- generate_profile_sequence(design$total_n, design$patterns, design$profile_prob)
  z <- absolute_minimization_assign(
    profile$X, pbc = design$pbc, weights = design$weights, seed = seed_allocation
  )
  A <- as.integer((z + 1L) / 2L)
  set.seed(normalise_seed(seed_outcome))
  n <- design$total_n
  if (model$type == "continuous_symmetric") {
    y <- rnorm(n, sd = model$shared_noise_sd) + 0.5 * z * model$delta[profile$id]
  } else if (model$type == "continuous_realistic") {
    y0 <- model$mu0[profile$id] + rnorm(n, sd = model$outcome_sd)
    eta <- rnorm(n, sd = model$individual_effect_sd)
    y <- y0 + A * (model$delta[profile$id] + eta)
  } else {
    probability <- ifelse(A == 1L, model$p1[profile$id], model$p0[profile$id])
    y <- rbinom(n, 1L, probability)
  }
  list(X = profile$X, profile_id = profile$id, z = z, A = A, y = y)
}

make_swift_direct_design <- function() {
  factor_prob <- c(0.466, 0.592, 0.287, 0.154, 0.300)
  patterns <- all_binary_patterns(5L)
  list(
    K = 5L,
    J = nrow(patterns),
    total_n = 404L,
    n_per_group = 202L,
    pbc = 0.80,
    weights = rep(1, 6L),
    patterns = patterns,
    factor_prob = factor_prob,
    profile_prob = profile_prob_independent(factor_prob, patterns),
    boundary = -0.12,
    alpha = 0.05
  )
}

make_swift_direct_model <- function(design) {
  coefficients <- c(-0.55, -0.25, -0.40, -0.25, -0.60)
  alpha0 <- calibrate_logistic_intercept(
    design$patterns, design$profile_prob, coefficients,
    interaction = 0, target_risk = 0.622
  )
  eta <- alpha0 + drop(design$patterns %*% coefficients)
  probability <- plogis(eta)
  list(
    type = "binary_swift",
    alpha0 = alpha0,
    coefficients = coefficients,
    p0 = probability,
    p1 = probability,
    effect = 0,
    d_at = function(boundary) probability - probability - boundary
  )
}

generate_swift_trial <- function(design, model, seed_profile, seed_allocation, seed_outcome) {
  set.seed(normalise_seed(seed_profile))
  profile <- generate_profile_sequence(design$total_n, design$patterns, design$profile_prob)
  z <- absolute_minimization_assign(
    profile$X, pbc = design$pbc, weights = design$weights, seed = seed_allocation
  )
  A <- as.integer((z + 1L) / 2L)
  set.seed(normalise_seed(seed_outcome))
  y <- rbinom(design$total_n, 1L, model$p0[profile$id])
  list(X = profile$X, profile_id = profile$id, z = z, A = A, y = y)
}
