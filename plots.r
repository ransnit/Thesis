source("forms.r")

plot.best.response <- function(p, ro)
{
  id <- function(x){return(x)}
  plot.function(id, from = 0, to = 1, lty = "dashed", ylab = "Best Response", xlab = "P_e")
  
  p <- c(p, 1)
  cond1 <- (p > ((ro - 1) / ro)) # True if x > ro
  
  if (cond1)
  {  
    segments(x0 = 0, y0 = 0, x1 = p[1], y1 = 0)
    segments(x0 = p[1], y0 = 1, x1 = p[2], y1 = 1)
    segments(x0 = p[2], y0 = 0, x1 = 1, y1 = 0)  
  }
  else
  {  
    segments(x0 = 0, y0 = 1, x1 = p[1], y1 = 1)
    segments(x0 = p[1], y0 = 0, x1 = p[2], y1 = 0)
    segments(x0 = p[2], y0 = 1, x1 = 1, y1 = 1)  
  }
  
  for(p0 in head(p, -1))
      segments(x0 = p0, y0 = 0, x1 = p0, y1 = 1, lty="dotted")
}

calc.and.plot <- function(gamma, nu, ro)
{
  cat("Calculating with parameters: gamma =", gamma,"; nu =", nu, "; ro =", ro, "\n")
  p <- calculate.pe(gamma, nu, ro)
  
  cat("Result: p_eq =", p, "\n")
  plot.best.response(p, ro)
}