is.prob <- function(p) { return ((p <=1) & (p >= 0)) }

find.roots <- function(a, b, c)
{
  stopifnot(a != 0)
#   if (a == 0)
#     return (-c / b)
    
  discriminant = b^2 - 4*a*c

  cat("Delta =", discriminant, "\n")

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
  
  cat("ro*(nu-1) =", ro*(nu-1), "and gamma+nu =", gamma+nu, "\n")
  
  x <- find.roots(a, b, c)
  cat("nu+ro*a-2*a =", nu+ro*a-2*a, "\n")
  cat("x_1, x_2 =", x, "\n")
  
  p <- (x - 1) / ro
  
  return (p)
}

calculate.pe2 <- function(gamma, ro)
{
  a = ro*(ro-1)
  b = ro*(1+gamma*(ro-2))
  c = gamma*(ro-1)
  
  if (a == 0)
    return (-c / b)
  
  return (find.roots(a,b,c))
}