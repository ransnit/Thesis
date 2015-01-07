MAX_RO <- 0.5*(1 + sqrt(5))
RO_RESO <- 0.05
RO_SEQ <- seq(RO_RESO, MAX_RO, RO_RESO)
ACC <- 1e-6
MIN_N <- 20
LP_VALS <- readRDS("LPVals.rds")

effective_ro <- function(ro, p) { return (ro*(1 - p/(1+p*ro)))}
fixedvals <- function(x) { return (x[!is.infinite(x) & !is.na(x)]) }

get.p.Lp <- function(ro)
{
  p <- LP_VALS[[paste(ro, "p")]]
  Lp <- LP_VALS[[paste(ro, "Lp")]]
  
  return (data.frame(p, Lp))
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

find.proper.range <- function(d)
{
  f <- function(n) { return (sum(head(d,n))) }
  d[d<0] <- 0
  d[d>1] <- 1
  s1 <- sapply(1:nrow(d), f)
  d <- d[which(s1 <= 1 & s1 >= 0), ]
  return(d)
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

plot.queue.length <- function(ro, accuracy = ACC, reso = 0.02)
{
  f <- function(p) { return (expected.queue.length(ro, p, accuracy)) }
  
  p <- seq(0, 1, reso)
  Lp <- sapply(p, f)
  
  plot(p, Lp, type='l')
  abline(v = ((ro - 1) / (ro * (2-ro))))
  
  return (cbind(p, Lp))
}

calc.p.eq <- function(ro, gamma, accuracy = ACC, reso = 0.02, pl = TRUE, readlp = TRUE)
{
  if (ro < 1 & gamma <= 1 - ro)
    return (0)
  
  f <- function(p) { return (expected.queue.length(ro, p, accuracy)) }
  
  if (readlp){
    d <- get.p.Lp(ro)
    p <- d$p
    Lp <- d$Lp
  }
  
  else{
    p <- seq(0, 1, reso)
    Lp <- sapply(p, f)
  }
  
  Vp <- gamma*Lp/(1+p*ro) - 1
  min_elem <- which.min(abs(Vp))
  
  if(pl)
  {
    Vs_p <- (1 + gamma*Lp*p*ro/(1+p*ro))
    Vn_p <- gamma*Lp
    ylimUp <- max(fixedvals(c(Vs_p, Vn_p, Vp)))
    ylimLo <- min(fixedvals(c(Vs_p, Vn_p, Vp)))
    plot(p, Vs_p, col="blue", type='l', xlab="", ylab="", xlim=c(0,1), ylim=c(ylimLo, ylimUp))
    textbox(x = c(0, 1), y = tail(fixedvals(Vs_p),1), textlist = c("Cost of Sensing"), box = F, justify = 'r' , cex = 0.7, col = "blue")
    par(new=T)
    plot(p, Vn_p, col="green", type='l', xlab="", ylab="", xlim=c(0,1), ylim=c(ylimLo, ylimUp))
    textbox(x = c(0, 1), y = tail(fixedvals(Vn_p),1), textlist = c("Cost of Not-Sensing"), box = F, justify = 'r' , cex = 0.7, col = "green")
    par(new=T)
    plot(p, Vp, col = ifelse( Vp == Vp[min_elem],'red','black'), type='b', lwd = 2, xlab="", ylab="", xlim=c(0,1), ylim=c(ylimLo, ylimUp))
    abline(0, 0, lty=2)
    title(paste("Value vs. p; gamma = ", gamma, "; ro = ", ro, sep = ""))
  }
  
  if (all(Vp[!is.na(Vp)] > 0))
    return (1)
  
  if (all(Vp[!is.na(Vp)] < 0))
    return (0)
  
  if (!is.na(p[min_elem]) & !is.infinite(p[min_elem]))
    return (p[min_elem])
  
  else
    return (NA)
}

p.ro.chart.for.given.gamma <- function(gamma, ro_seq = RO_SEQ)
{
  calc.p.eq.wrapper <- function(ro) { return (calc.p.eq(ro, gamma, reso = 0.01, pl = FALSE)) }
  p_vals <- sapply(ro_seq, calc.p.eq.wrapper)
  plot(ro_seq, p_vals, type='b')
  return (p_vals)
}

require('plotrix')

accumulated.p.ro.chart <- function(p_vec_list, gamma_vals, ro_seq = RO_SEQ, ro_xlim = c(0, MAX_RO))
{
  ngraphs <- length(gamma_vals)
  colors <- c("blue", "red", "green")
  for (i in 1:ngraphs)
  {
    p_vals <- p_vec_list[[i]]
    gamma <- gamma_vals[i]
    color <- colors[(i%%length(colors))+1]
    plot(ro_seq, p_vals, xlab = "", ylab = "", xlim = ro_xlim, ylim = c(0,1), type = 'b', lwd = 1, pch = 20, col = color)
    txtanchor <- which(p_vals >= 0.1*ngraphs-0.075*(ngraphs-i))[1]
    tx <- ro_seq[txtanchor]
    ty <- p_vals[txtanchor]
    textbox(x = c(tx-1, tx+1), y = ty+0.025, textlist = c("1/gamma =",1/gamma), box = F, justify = 'c' , cex = 0.7, col = color)
    par(new = T)
  }
  title(main = "P(eq) vs. ro", xlab = "ro", ylab = "P(eq)")
  par(new = F)
}

p.ro.chart <- function(vals, ro_xlim = c(0, MAX_RO))
{
  p_vals_list <- readRDS("p_vals.rds")
  gamma_vals <- readRDS("gamma_vals.rds")
  indices <- match(vals, gamma_vals)
  
  stopifnot(all(!is.na(indices)))
  
  accumulated.p.ro.chart(p_vals_list[indices], gamma_vals[indices], RO_SEQ, ro_xlim)
}  