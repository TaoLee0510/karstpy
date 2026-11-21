# Generate reference fits from the isolated R karst 1.1.3 installation.
args <- commandArgs(TRUE)
native_only <- length(args) == 3L && identical(args[[3L]], "--native-only")
if (length(args) != 2L && !native_only) {
  stop("Usage: Rscript tools/parity/generate.R R_LIBRARY OUTPUT_DIR [--native-only]")
}
library(karst, lib.loc = normalizePath(args[[1L]]))
stopifnot(packageVersion("karst") == "1.1.3", packageVersion("TMB") == "1.9.19")
dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[[2L]])
cases <- list()

save_native <- function(name, fit) {
  p <- karst:::sparse_joint_problem(fit$data, fit$settings)
  joint <- TMB::MakeADFun(p$tmb_data, p$parameters, DLL = "karst", silent = TRUE)
  marginal <- TMB::MakeADFun(p$tmb_data, p$parameters, map = p$map, random = p$random,
                            DLL = "karst", silent = TRUE,
                            inner.control = list(maxit = 200, tol = 1e-10))
  set.seed(3L)
  points <- lapply(c(0, 0.002, 0.01), function(sd) joint$par + rnorm(length(joint$par), 0, sd))
  tests <- lapply(points, function(point) list(point = point, fn = joint$fn(point),
                                               gr = joint$gr(point), he = joint$he(point)))
  laplace <- lapply(c(0, 0.01, -0.015), function(step) {
    point <- marginal$par + step
    value <- marginal$fn(point)
    best <- marginal$env$last.par.best
    mode <- marginal$env$last.par
    gradient <- marginal$gr(point)
    list(fixed = point, value = value, grad = gradient, best = unname(best),
         mode = unname(mode), gradient_mode = unname(marginal$env$last.par))
  })
  fixture <- list(model = p$tmb_data$model, data = p$tmb_data[setdiff(names(p$tmb_data), "model")],
                  parameters = p$parameters, random = p$random, free = p$free_scale,
                  member = p$tmb_data$member, joint = tests, laplace = laplace)
  jsonlite::write_json(fixture, file.path(out, paste0("native-", name, ".json")),
                       # NA uses 15 significant digits in jsonlite 2.0.0.
                       digits = 17L, auto_unbox = TRUE)
}

save_fit <- function(name, case, seed, fit, counts, times, history, series, beta, birth_rate = 1, rate_group = NULL) {
  utils::write.csv(counts, file.path(out, paste0(name, ".counts.csv")))
  utils::write.csv(fit$summary, file.path(out, paste0(name, ".summary.csv")))
  d <- fit$diagnostics
  cases[[length(cases) + 1L]] <<- list(
    name = name, case = case, seed = seed, times = times, history = history,
    series = series, beta = beta, birth_rate = birth_rate, rate_group = rate_group,
    selected_history = d$history, convergence = d$convergence, objective = d$objective,
    posterior_objective = d$posterior_objective, dropped = fit$dropped, left_out = d$dropped_samples$sample
  )
}

for (seed in seq_len(if (native_only) 1L else 10L)) {
  landscape <- simulate_grf_landscape(n_chr = 3, min_cn = 1, max_cn = 4, seed = seed)
  short <- simulate_sparse_counts(landscape, beta = 0.001, dt = 0.3, seed = 200L + seed,
                                 n0 = 220, n1 = 320, initial_population = 600L)$counts
  y <- simulate_sparse_counts(landscape, beta = 0.001, dt = 2, seed = 400L + seed,
                             n0 = 220, n1 = 320, initial_population = 600L)$counts
  set.seed(10L + seed)
  middle <- vapply(c(1/3, 2/3), function(w) {
    as.integer(rmultinom(1L, 280L, (1 - w) * y[, 1L] / sum(y[, 1L]) + w * y[, 2L] / sum(y[, 2L])))
  }, integer(nrow(y)))
  series_counts <- cbind(y[, 1L], middle, y[, 2L])
  specs <- list(pair = list(y, c(0, 2), c(2, 5, 10), "anchored"),
                short = list(short, c(0, 0.3), c(0.3, 2, 8), "anchored"),
                anchored = list(series_counts, c(0, .65, 1.35, 2), c(2, 5, 10), "anchored"),
                pairwise = list(series_counts, c(0, .65, 1.35, 2), c(2, 5, 10), "pairwise"),
                state_space = list(series_counts, c(0, .65, 1.35, 2), c(2, 5, 10), "state_space"))
  for (case in names(specs)) {
    if (native_only && case == "short") next
    spec <- specs[[case]]
    fit <- suppressWarnings(fit_sparse_joint(spec[[1L]], times = spec[[2L]], beta = 0.001, max_cn = 4,
                                             history = spec[[3L]], series = spec[[4L]]))
    name <- paste(case, seed, sep = "-")
    if (!native_only) save_fit(name, case, seed, fit, spec[[1L]], spec[[2L]], spec[[3L]], spec[[4L]], 0.001)
    if (seed == 1L && case != "short") save_native(case, fit)
    cat(name, "\n")
  }
  if (native_only) next
  rates <- suppressWarnings(fit_sparse_joint_rates(y, dt = 2, beta = 0.001, history = 2, max_cn = 4,
                                                   birth_rate_range = c(0.6, 1.2)))
  for (i in seq_along(rates$rates)) {
    name <- paste("rates", seed, i, sep = "-")
    save_fit(name, "rates", seed, rates$fits[[i]], y, c(0, 2), 2, "anchored", 0.001,
             birth_rate = rates$rates[[i]], rate_group = i)
    cat(name, "\n")
  }
  sim <- simulate_sparse_counts(landscape, beta = 0.01, dt = 2, seed = seed, n0 = 90, n1 = 110,
                                 initial_population = 250L, overdispersion = if (seed %% 2) 0 else 70)
  write.csv(sim$counts, file.path(out, paste0("sim-", seed, ".csv")))
  writeBin(sim$truth_observed, file.path(out, paste0("sim-", seed, ".bin")), endian = "little")
}
if (!native_only) jsonlite::write_json(cases, file.path(out, "cases.json"), digits = 17L, auto_unbox = TRUE)
