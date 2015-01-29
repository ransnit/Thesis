source("probs.r")

gamma_vals <- unique(c(seq(1, 100, 1), seq(1, 100, 1)^(-1)))
# gamma_vals <- c(1, 2)
p_vals <- lapply(gamma_vals, p.ro.chart.for.given.gamma)

saveRDS(gamma_vals, "gamma_vals.rds")
saveRDS(p_vals, "p_vals.rds")

qls <- function(ro)
{
  f <- function(p) { return (expected.queue.length(ro, p, ACC)) }
  f0 <- function(p) { return (conditional.expected.queue.length(ro, p, 0, ACC)) }
  p <- seq(0, 1, 0.01)
  Lp <- sapply(p, f)
  Lp0 <- sapply(p, f0)
  return (data.frame(p, Lp, Lp0))
}

RO_SEQ <- c(0.4, 0.5, 0.6) #TODO!!!

l <- sapply(RO_SEQ, qls)
names(l) <- paste(RO_SEQ)
saveRDS(l, "LPVals.rds")

