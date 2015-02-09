MAX_RO <- 0.5*(1 + sqrt(5))
RO_RESO <- 0.05
RO_SEQ <- seq(RO_RESO, MAX_RO, RO_RESO)
ACC <- 1e-6
MIN_N <- 20
LP_VALS <- readRDS("LPVals.rds")
CEX_MAIN <- 2.5
CEX_LAB <- 1.4

effective_ro <- function(ro, p) { return (ro*(1 - p/(1+p*ro)))}
fixedvals <- function(x) { return (x[!is.infinite(x) & !is.na(x)]) }

get.p.Lp <- function(ro)
{
  p <- LP_VALS[[paste(ro, "p")]]
  Lp <- LP_VALS[[paste(ro, "Lp")]]
  Lp0 <- LP_VALS[[paste(ro, "Lp0")]]
  
  return (data.frame(p, Lp, Lp0))
}

calc.prob.table.by.p00 <- function(ro, p, accuracy = ACC)
{
  n <- max(MIN_N, ceiling(log(accuracy, base = effective_ro(ro, p))))
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

calc.p00.from.table <- function(d, accuracy)
{
  a <- sum(d$P_n0_a, d$P_n1_a)
  b <- sum(d$P_n0_b, d$P_n1_b)
  
  return ((1 - b) / a)
}

calc.p00 <- function(ro, p, accuracy = ACC)
{
  d <- calc.prob.table.by.p00(ro, p, accuracy)
  return (calc.p00.from.table(d, accuracy))
}

calc.prob.table <- function(d, accuracy)
{
  p00 <- calc.p00.from.table(d, accuracy)
  
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

calc.conditional.expected.value <- function(d, k)
{
  n <- nrow(d)
  p <- d[, paste("P_n", k, sep = "")]
  pn <- p / sum(p)
  t <- (pn * seq(0,n-1))
  return (sum(t))
}

find.proper.range <- function(d)
{
  f <- function(n) { return (sum(head(d,n))) }
  d[d<0] <- 0
  d[d>1] <- 1
  s1 <- sapply(1:nrow(d), f)
  d <- d[which(s1 <= 1 & s1 >= 0), ]
  return(d)
}

prob.n <- function(n, ro, p, accuracy = ACC)
{
  if (ro >= MAX_RO)
    return (0)
  
  if (p <=  ((ro - 1) / (ro * (2-ro))))
    return (0)
  
  if (p == 0)
    return ((1-ro)*ro^n)
  
  d <- calc.prob.table(calc.prob.table.by.p00(ro, p, accuracy), accuracy)
  d <- find.proper.range(d)
  
  if(nrow(d) < 1)
    return (NA)
  
  return (d$P_n0[n+1]+d$P_n1[n+1])
}

prob.n0 <- function(n, ro, p, accuracy = ACC)
{
  if (ro >= MAX_RO)
    return (0)
  
  if (p <=  ((ro - 1) / (ro * (2-ro))))
    return (0)
  
  if (p == 0)
    return ((1-ro)*ro^n)
  
  d <- calc.prob.table(calc.prob.table.by.p00(ro, p, accuracy), accuracy)
  d <- find.proper.range(d)
  
  if(nrow(d) < 1)
    return (NA)
  
  return (d$P_n0[n+1])
}

probs.to.n <- function(n, ro, p, accuracy = ACC)
{
  if (ro >= MAX_RO)
    return (NA)
  
  if (p <=  ((ro - 1) / (ro * (2-ro))))
    return (NA)
  
  if (p == 0)
    return ((1-ro)*(ro))
  
  d <- calc.prob.table(calc.prob.table.by.p00(ro, p, accuracy), accuracy)
  d <- find.proper.range(d)
  
  if(nrow(d) < 1)
    return (NA)
  
  return (d$P_n0[1]+d$P_n1[1]- p*sum(d$P_n0[1:(n+1)]))
}

expected.queue.length <- function(ro, p, accuracy = ACC)
{
  if (ro >= MAX_RO)
    return (Inf)
  
  if (p <=  ((ro - 1) / (ro * (2-ro))))
    return (Inf)
  
  if (p == 0)
    return (ro / (1-ro))
  
  d <- calc.prob.table(calc.prob.table.by.p00(ro, p, accuracy), accuracy)
  d <- find.proper.range(d)
  
  if(nrow(d) < 1)
    return (NA)
  
  residue <- 1 - sum(d)
  
  l <- calc.expected.value(d) + residue*nrow(d)
  mm1l <- (1-p)*ro / (1 - (1-p)*ro)
  
  if (l <= mm1l)
    return (NA)
  
  return (l)
}

conditional.expected.queue.length <- function(ro, p, k, accuracy = ACC)
{
  if (ro >= MAX_RO)
    return (Inf)
  
  if (p <=  ((ro - 1) / (ro * (2-ro))))
    return (Inf)
  
  if (p == 0)
    return (ro / (1-ro))
  
  d <- calc.prob.table(calc.prob.table.by.p00(ro, p, accuracy), accuracy)
  d <- find.proper.range(d)
  
  if(nrow(d) < 1)
    return (NA)
  
  residue <- ((1-k)+k*p*ro)/(p*ro+1) - sum(d[, paste("P_n", k, sep = "")])
  l <- calc.conditional.expected.value(d, k) + residue*nrow(d)
  
  return (l)
}

conditional.expected.queue.length2 <- function(ro, p, accuracy = ACC)
{
  p_00 <- calc.p00(ro,p)
  g_ <- 1 - ro*(1 - p / (1 + p*ro))
  dg_ <- g_ + (1-ro)*p_00
  g <- -p*ro^2 -ro +2*p*ro +1
  dg <- ro^2 -3*p*ro^2 -4*ro + 4*p*ro +2
  res <- (dg_*g - dg*g_) / g^2
  
  return ((1+p*ro)*res)
}

plot.queue.length <- function(ro, accuracy = ACC, reso = 0.02)
{
  f1 <- function(p) { return (expected.queue.length(ro, p, accuracy)) }
  f2 <- function(p) { return (conditional.expected.queue.length(ro, p, 1, accuracy)) }
  f3 <- function(p) { return (conditional.expected.queue.length(ro, p, 0, accuracy)) }
  
  p <- seq(0, 1, reso)
  Lp1 <- sapply(p, f1)
  Lp2 <- sapply(p, f2)
  Lp3 <- sapply(p, f3)
  
  ylim <- c(min(fixedvals(c(Lp1, Lp2, Lp3))), max(fixedvals(c(Lp1, Lp2, Lp3))))
  
  plot(p, Lp1, type='l', ylim = ylim, ylab = "")
  par(new=T)
  plot(p, Lp2, type='l', col = "blue", ylim = ylim, ylab = "")
  par(new=T)
  plot(p, Lp3, type='l', col = "red", ylim = ylim, ylab = "")
  abline(v = ((ro - 1) / (ro * (2-ro))))
  
#   return (cbind(p, Lp))
}

plot.prob.n <- function(n, ro, accuracy = ACC, reso = 0.02)
{
  f <- function(p) { return (prob.n(n, ro, p, accuracy)) }
  
  p <- seq(0, 1, reso)
  pn <- sapply(p, f)
  
  plot(p, pn, type='l')
  
  return (pn)
}

plot.prob.n0 <- function(n, ro, accuracy = ACC, reso = 0.02)
{
  f <- function(p) { return (prob.n0(n, ro, p, accuracy)) }
  
  p <- seq(0, 1, reso)
  pn <- sapply(p, f)
  
  plot(p, pn, type='l')
  
  return (pn)
}

plot.probs.to.n <- function(n, ro, accuracy = ACC, reso = 0.02)
{
  f <- function(p) { return (probs.to.n(n, ro, p, accuracy)) }
  
  p <- seq(0, 1, reso)
  pn <- sapply(p, f)
  
  plot(p, pn, type='l')
  
  return (pn)
}

calc.p.eq <- function(ro, gamma, accuracy = ACC, reso = 0.02, readlp = TRUE)
{
  if (ro <= 1 / (1 + gamma))
    return (0)
  
  if (readlp){
    d <- get.p.Lp(ro)
    p <- d$p
    Lp <- d$Lp
    Lp0 <- d$Lp0
  }
  
  else{
    f <- function(p) { return (expected.queue.length(ro, p, accuracy)) }
    f0 <- function(p) { return (conditional.expected.queue.length(ro, p, 0, accuracy)) }
    
    p <- seq(0, 1, reso)
    Lp <- sapply(p, f)
    Lp0 <- sapply(p, f0)
  }
  
  Vp <- gamma*Lp0/(1+p*ro) - 1
  min_elem <- which.min(abs(Vp))
  
  if (all(Vp[!is.na(Vp)] < 0))
    return (0)
  
  if (all(Vp[!is.na(Vp)] > 0))
    return (1)
  
  if (!is.na(p[min_elem]) & !is.infinite(p[min_elem]))
    return (p[min_elem])
  
  else
    return (NA)
}

plot.costs.chart <- function(ro, gamma, accuracy = ACC, reso = 0.02, readlp = TRUE)
{
  if (readlp){
    d <- get.p.Lp(ro)
    p <- d$p
    Lp <- d$Lp
    Lp0 <- d$Lp0
  }
  
  else{
    f <- function(p) { return (expected.queue.length(ro, p, accuracy)) }
    f0 <- function(p) { return (conditional.expected.queue.length(ro, p, 0, accuracy)) }
    
    p <- seq(0, 1, reso)
    Lp <- sapply(p, f)
    Lp0 <- sapply(p, f0)
  }
  
  Vs_p <- (1 + gamma*Lp0*p*ro/(1+p*ro))
  Vn_p <- gamma*Lp
  ylimUp <- max(fixedvals(c(Vs_p, Vn_p)))
  ylimLo <- min(fixedvals(c(Vs_p, Vn_p)))
  plot(p, Vs_p, col="gray70", type='l', xlab="", ylab="", xlim=c(0,1), ylim=c(ylimLo, ylimUp), lwd = 2)
  text(x = 0.95, y = tail(fixedvals(Vs_p),1)+0.1, labels = bquote(C[S]), cex = 1.7, col = "gray70")
  par(new=T)
  plot(p, Vn_p, col="gray10", type='l', xlab="", ylab="", xlim=c(0,1), ylim=c(ylimLo, ylimUp), lwd = 2)
  text(x = 0.95, y = tail(fixedvals(Vn_p),1)+0.1, labels = bquote(C[N]), cex = 1.7, col = "gray10")
  par(new=T)
  title(main = bquote("Cost v.s." ~ p ~ ";" ~ gamma == ~.(gamma) ~ "," ~ rho == ~ .(ro)), cex.main = CEX_MAIN,
        ylab = "Cost", cex.lab = CEX_LAB, xlab = bquote(p))
}

p.ro.chart.for.given.gamma <- function(gamma, ro_seq = RO_SEQ)
{
  calc.p.eq.wrapper <- function(ro) { return (calc.p.eq(ro, gamma, reso = 0.01)) }
  p_vals <- sapply(ro_seq, calc.p.eq.wrapper)
#   plot(ro_seq, p_vals, type='b')
  return (p_vals)
}

require('plotrix')

accumulated.p.ro.chart <- function(p_vec_list, gamma_vals, ro_seq = RO_SEQ, ro_xlim = c(0, MAX_RO))
{
  ngraphs <- length(gamma_vals)
  colors <- c("gray0", "gray45", "gray70")
  for (i in 1:ngraphs)
  {
    p_vals <- p_vec_list[[i]]
    gamma <- gamma_vals[i]
    color <- colors[(i%%length(colors))+1]
    plot(ro_seq, p_vals, xlab = "", ylab = "", xlim = ro_xlim, ylim = c(0,1), type = 'b', lwd = 1, pch = 20, col = color)
    txtanchor <- which(p_vals >= 0.1*ngraphs-0.075*(ngraphs-i))[1]
    tx <- ro_seq[txtanchor]
    ty <- p_vals[txtanchor]
#     textbox(x = c(tx-1, tx+1), y = ty+0.025, textlist = c("1/gamma =",round(1/gamma,4)), box = F, justify = 'c' , cex = 0.7, col = color)
    tmp <- round(1/gamma,4)
    text(x = tx, y = ty+0.025, labels = bquote(gamma^-1  == ~ .(tmp)), cex = 1.2, col = color)
    par(new = T)
  }
  title(main = bquote(p[e] ~ v.s. ~ rho), cex.main = CEX_MAIN, 
        xlab = bquote(rho), ylab = bquote(p[e]), cex.lab = CEX_LAB)
  par(new = F)
}

temp <- function(ro)
{
  l1<- sapply(seq(0, 1, 0.05), FUN = function(p){return(conditional.expected.queue.length(ro,p,0))})
  l2<- sapply(seq(0, 1, 0.05), FUN = function(p){return(conditional.expected.queue.length2(ro,p))})
  
  ylim <- c(0, max(fixedvals(c(l1,l2))))
  
  plot(l1, type="l", col="red", ylim = ylim)
  par(new=T)
  plot(l2, type="l", col="blue", ylim = ylim)
}

p.ro.chart <- function(vals, ro_xlim = c(0, MAX_RO))
{
  p_vals_list <- readRDS("p_vals.rds")
  gamma_vals <- readRDS("gamma_vals.rds")
  indices <- match(vals, gamma_vals)
  
  stopifnot(all(!is.na(indices)))
  
  accumulated.p.ro.chart(p_vals_list[indices], gamma_vals[indices], RO_SEQ, ro_xlim)
}  