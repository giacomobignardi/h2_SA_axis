# Profile-likelihood confidence interval for a genetic correlation in the twin models (R2.5)
# rA is fixed at a value r0 (all other parameters re-estimated, starting from the unconstrained
# solution); each bound is the r0 at which the likelihood-ratio statistic (delta chi-square, 1 df)
# reaches the critical value, located by root finding. If the statistic stays below the critical
# value up to |r0| = edge, the interval is reported as reaching the edge of the admissible range.
profile_ci = function(model, data, cov_lab, var1, var2, r_lab,
                      level = .95, tol = .001, edge = .99) {
  fit_m = function(m, ...) sem(m, data = data, group = "zyg", group.label = c(1, 2),
                               missing = "ML", std.ov = T, ...)
  base = fit_m(model)
  wald = parameterEstimates(base, level = level) %>% filter(label == r_lab)
  chi0 = fitMeasures(base, "chisq")
  crit = qchisq(level, df = 1)                                   # 3.84 for a 95% interval
  dchi = function(r0) {                                          # delta chi-square with rA fixed at r0
    f = tryCatch(suppressWarnings(fit_m(paste0(model, sprintf(
      "\n %s == %.6f*sqrt(%s*%s)\n", cov_lab, r0, var1, var2)), start = base)),
      error = function(e) NULL)
    if (is.null(f) || !lavInspect(f, "converged")) NA_real_ else unname(fitMeasures(f, "chisq") - chi0)
  }
  find_bound = function(est, to) {                               # search between estimate and edge
    d_edge = dchi(to)
    if (is.na(d_edge)) return(list(bound = NA_real_, note = "refit failed at edge"))
    if (d_edge <= crit) return(list(bound = to, note = "reaches edge"))
    r = tryCatch(uniroot(function(x) dchi(x) - crit, sort(c(est, to)), tol = tol)$root,
                 error = function(e) NA_real_)
    list(bound = r, note = if (is.na(r)) "search failed" else "")
  }
  est = wald$est
  if (abs(est) >= edge) {
    lo = hi = list(bound = NA_real_, note = "estimate outside admissible range")
  } else {
    lo = find_bound(est, -edge); hi = find_bound(est, edge)
  }
  tibble(r             = r_lab,
         est           = est,
         wald_lower    = wald$ci.lower, wald_upper = wald$ci.upper,
         profile_lower = lo$bound,      profile_upper = hi$bound,
         note_lower    = lo$note,       note_upper    = hi$note)
}