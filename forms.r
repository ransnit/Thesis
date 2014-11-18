is.prob <- function(p) { return ((p <=1) & (p >= 0)) }

find.roots <- function(a, b, c)
{
  stopifnot(a != 0)
#   if (a == 0)
#     return (-c / b)
    
  discriminant = b^2 - 4*a*c

  stopifnot(discriminant >= 0)
  
  x1 = (-b - sqrt(discriminant)) / (2*a)
  x2 = (-b + sqrt(discriminant)) / (2*a)
  
  return (c(x1, x2))
}

calculate.pe <- function(gamma, nu, ro)
{
  a = gamma + nu
  b = -(nu + ro*gamma + nu*ro)
  c = nu*ro - ro
  
  p <- (find.roots(a, b, c) - 1) / ro
  
#  return (p[is.prob(p)])
  return (p)
}
