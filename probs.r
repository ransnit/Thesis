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
  t <- (d$P_n0+d$P_n1) * seq(0,n-1)
  return (sum(t))
}

find.proper.range <- function(d)
{
  f <- function(n) { return (sum(head(d,n))) }
  is.prob <- function(v) { return (all(v >= 0 & v <= 1)) }
  s1 <- sapply(1:nrow(d), f)
  d <- d[which(s1 <= 1 & s1 > 0), ]
  s2 <- apply(X=d, MARGIN = 1, FUN = is.prob)
  return (d[s2,])
}

expected.queue.length <- function(ro, p, n)
{
  d <- calc.prob.table(calc.prob.table.by.p00(ro, p, n))
  d <- find.proper.range(d)
  return (calc.expected.value(d))
}

prob.as.a.function.of.p00.and.p01 <- function(p, ro, n)
{
  sumlist <- function(l) { return(Reduce("+",l)) }
  lp_n0 <- list(c(1,0))
  lp_n1 <- list(c(0,1))
  
  if (n > 0)
  {
    for(i in 2:(n+1))
    {
      s <- sumlist(head(lp_n1, i-1)) - ro*p*sumlist(head(lp_n0, i-1))
      lp_n1[[i]] <- ro*lp_n1[[i-1]] + s
      lp_n0[[i]] <- (1-p)*ro*lp_n0[[i-1]] - s
    }
  }
  
  cat("P_",n,"0 = ", sep="")
  show(lp_n0[[n+1]])
  cat("P_",n,"1 = ", sep="")
  show(lp_n1[[n+1]])
  
  return (lp_n0[[n+1]]+lp_n1[[n+1]])
}

p11.and.p10 <- function(p, ro)
{
  return ( ro*c(1-p, 1) )
}

p21.and.p20 <- function(p, ro)
{
  t1 <- ro - 2*p*ro
  t2 <- p + ro
  return ( ro*c(t1, t2) )
}

p31.and.p30 <- function(p, ro)
{
  t1 <- ro^2 - 3*p*ro^2 - 2*p^2*ro
  t2 <- 3*p*ro + 2*p + ro^2
  return ( ro*c(t1, t2) )
}

p41.and.p40 <- function(p, ro)
{
  t1 <- ro^3 - 4*p*ro^3 - 7*p^2*ro^2 - 4*p^2*ro
  t2 <- 7*p*ro + 2*p^2*ro + 6*p*ro^2 + 4*p + ro^3
  return ( ro*c(t1, t2) )
}

p51.and.p50 <- function(p, ro)
{
  t1 <- ro^4 - 16*p^2*ro^2 - 5*p*ro^4 - 16*p^2*ro^3 - 4*p^3*ro^2 - 8*p^2*ro
  t2 <- 16*p*ro^2 + 16*p*ro + ro^4 + 10*p*ro^3 + 9*p^2*ro^2 + 8*p^2*ro + 8*p
  return ( ro*c(t1, t2) )
}

MAX_RO <- 0.5*(1 + sqrt(5))
RO_RESO <- 0.05
RO_SEQ <- seq(RO_RESO, MAX_RO, RO_RESO)

calc.p.eq <- function(ro, gamma, n = 30, reso = 0.005, pl = TRUE)
{
  if (ro >= MAX_RO)
    return (NA)
  
  if (ro < 1 & gamma <= 1 - ro)
    return (0)
  
  f <- function(p)
  {
    if ((1-p)*ro + (p*ro)^2/(1+p*ro) >= 1)
      return (Inf)
    Lp <- expected.queue.length(ro, p, n)
    return (gamma*(Lp - 1) - p*ro - 1)
  }
  
  p <- seq(reso, 1, reso)
  V_p <- sapply(p, f)
  min_elem <- which.min(abs(V_p))
  
  if(pl)
  {
    plot(p, V_p, col = ifelse( V_p == V_p[min_elem],'red','black'), lwd = 2)
    abline(0, 0, lty=2)
    title(paste("Value vs. p; gamma = ", gamma, "; ro = ", ro, sep = ""))
  }
  
  return (p[min_elem])
}

p.ro.chart.for.given.gamma <- function(gamma, ro_seq = RO_SEQ)
{
  calc.p.eq.wrapper <- function(ro) { return (calc.p.eq(ro, gamma, reso = 0.01, pl = FALSE)) }
  p_vals <- sapply(ro_seq, calc.p.eq.wrapper)
  plot(ro_seq, p_vals, type='b')
  return (p_vals)
}

require('plotrix')

accumulated.p.ro.chart <- function(p_vec_list, gamma_vals, ro_seq = RO_SEQ)
{
  ngraphs <- length(gamma_vals)
  colors <- c("blue", "red", "green")
  for (i in 1:ngraphs)
  {
    p_vals <- p_vec_list[[i]]
    gamma <- gamma_vals[i]
    color <- colors[(i%%length(colors))+1]
    plot(ro_seq, p_vals, xlab = "", ylab = "", xlim = c(0, MAX_RO), ylim = c(0,1), type = 'b', lwd = 1, pch = 20, col = color)
    txtanchor <- which(p_vals >= 0.1*ngraphs-0.075*(ngraphs-i))[1]
    tx <- ro_seq[txtanchor]
    ty <- p_vals[txtanchor]
    textbox(x = c(tx-1, tx+1), y = ty+0.025, textlist = c("1/gamma =",1/gamma), box = F, justify = 'c' , cex = 0.7, col = color)
    par(new = T)
  }
  title(main = "P(eq) vs. ro", xlab = "ro", ylab = "P(eq)")
  par(new = F)
}

p.ro.chart <- function(vals)
{
  p_vals_list <- readRDS("p_vals.rds")
  gamma_vals <- readRDS("gamma_vals.rds")
  indices <- match(vals, gamma_vals)
  
  stopifnot(all(!is.na(indices)))
  
  accumulated.p.ro.chart(p_vals_list[indices], gamma_vals[indices], RO_SEQ)
}  