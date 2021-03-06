source("forms.r")

bound.val.prob <- function(p)
{
  p <- max(p, 0)
  p <- min(p, 1)
  return (p)
}

seg.bounds <- function(v)
{
  stopifnot(length(v) > 1)
  
  res <- matrix(data = NA, nrow = (length(v)-1), ncol = 2)
  res[,1] <- head(v, -1)
  res[,2] <- tail(v, -1)
  return(res)
}

plot.best.response <- function(p, gamma, nu, ro)
{
  id <- function(x){return(x)}
  plot.function(id, from = 0, to = 1, lty = "dashed", ylab = "Best Response", xlab = "p_eq")
  
  tmp <- sapply(c(0, p, 1), bound.val.prob)
  tmp <- sort(unique(tmp))
  segs <- seg.bounds(tmp)
  
  draw.seg <- function(seg)
  {
    mid <- mean(seg)
    flip <- ((nu + gamma) > 0)
    br_noflip <- (p[1] < mid) && (mid < p[2])
    br <- as.integer(xor(flip, br_noflip))
    segments(x0 = seg[1], y0 = br, x1 = seg[2], y1 = br)
  }
  
  apply(X = segs, MARGIN = 1, FUN = draw.seg)
  
  for(p0 in p[is.prob(p)])
      segments(x0 = p0, y0 = 0, x1 = p0, y1 = 1, lty = "dotted")
}

calc.and.plot <- function(gamma, nu, ro)
{
  cat("Calculating with parameters: gamma =", gamma,"; nu =", nu, "; ro =", ro, "\n")
  p <- calculate.pe(gamma, nu, ro)
  
	cat("Roots are:", p, "\n")
	plot.best.response(p, gamma, nu, ro)
  
  cat("Result: p_eq =", p[is.prob(p)], "\n")
}

calc.and.plot.wrapper <- function(lambda, mu, c_w, c_s, R)
{
  gamma = (c_s*mu) / c_w
  nu = (R*mu) / c_w
  ro = lambda / mu
  
  calc.and.plot(gamma, nu, ro)
}
