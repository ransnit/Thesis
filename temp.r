source("probs.r")

gamma_vals <- unique(c(seq(1, 100, 1), seq(1, 100, 1)^(-1)))
# gamma_vals <- c(1, 2)
p_vals <- lapply(gamma_vals, p.ro.chart.for.given.gamma)

saveRDS(gamma_vals, "gamma_vals.rds")
saveRDS(p_vals, "p_vals.rds")

qls <- function(ro)
{
  f <- function(p) { return (expected.queue.length(ro, p, ACC)) }
  p <- seq(0, 1, 0.01)
  Lp <- sapply(p, f)
  return (data.frame(p, Lp))
}

# l <- sapply(RO_SEQ, qls)
# names(l) <- paste(RO_SEQ)
# saveRDS(l, "LPVals.rds")

