calc.prob.table.by.p00 <- function(ro, p, n)
{
  d <- data.frame(matrix(NA, n, 6))
  P_n0.l <- c("P_n0_a", "P_n0_b")
  P_n1.l <- c("P_n1_a", "P_n1_b")
  S_n.l <- c("S_n_a", "S_n_b")
  names(d) <- c(P_n0.l, P_n1.l, S_n.l)
  
  t <- 1 - ro*(1 - p / (1 + p*ro)) # t is p00+p01
  
  d[1,] <- c(1, 0, -1, t, NA, NA)
  d[1, S_n.l] <- d[1, P_n1.l] - ro*p*d[1, P_n0.l]
  
  for(i in 1:(n-1))
  {
    d[i+1, P_n1.l] <- ro*d[i, P_n1.l] + d[i, S_n.l]
    d[i+1, P_n0.l] <- (1-p)*ro*d[i, P_n0.l] - d[i, S_n.l]
    d[i+1, S_n.l] <- d[i, S_n.l] + d[i+1, P_n1.l] - ro*p*d[i+1, P_n0.l]
  }
  
  return (d)
}

calc.p00.from.table <- function(d)
{
  a <- sum(d$P_n0_a, d$P_n1_a)
  b <- sum(d$P_n0_b, d$P_n1_b)
  
  return ((1 - b) / a)
}

calc.prob.table <- function(d)
{
  p00 <- calc.p00.from.table(d)
  
  P_n0 <- d$P_n0_a*p00 + d$P_n0_b
  P_n1 <- d$P_n1_a*p00 + d$P_n1_b

  t <- data.frame(cbind(P_n0, P_n1))
  names(t) <- c("P_n0", "P_n1")
  
  return (t)
}

calc.expected.value <- function(d)
{
  n <- nrow(d)
  t <- d$P_n0+d$P_n1 * seq(0,n-1)
  return (sum(t))
}

find.proper.range <- function(d)
{
  f <- function(n) { return (sum(head(d,n))) }
  s <- sapply(1:nrow(d), f)
  return (d[which(s <= 1),])
}