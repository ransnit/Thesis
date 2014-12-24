source("probs.r")

gamma_vals <- c(1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32, 5, 10, 20, 1/5, 1/10, 1/20)
p_vals <- lapply(gamma_vals, p.ro.chart.for.given.gamma)

saveRDS(gamma_vals, "gamma_vals.rds")
saveRDS(p_vals, "p_vals.rds")