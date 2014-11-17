find.roots <- function(a, b, c)
{
  delta = b^2 - 4*a*c
  x1 = (-b + sqrt(delta)) / (2*a)
  x2 = (-b - sqrt(delta)) / (2*a)
  
  return (c(x1, x2))
}

calculate.pe <- function(gamma, nu, ro)
{
  a = gamma + nu
  b = -(nu + ro*gamma + nu*ro)
  c = nu*ro - ro
  
  p <- (find.roots(a, b, c) - 1) / ro
  
  is.prob <- function(p) { return ((p <=1) & (p >= 0)) }
  
  return (p[is.prob(p)])
}