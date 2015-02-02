source("probs.r")

gamma_vals <- unique(c(seq(1, 100, 1), seq(1, 100, 1)^(-1)))
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
# 
# l <- sapply(RO_SEQ, qls)
# 
# rsp <- paste(RO_SEQ, "p")
# rslp <- paste(RO_SEQ, "Lp")
# rslp0 <- paste(RO_SEQ, "Lp0")
# 
# names(l) <- c(t(cbind(rsp, rslp, rslp0)))
# saveRDS(l, "LPVals.rds")

