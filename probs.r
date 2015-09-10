Phi <- 0.5*(1 + sqrt(5))
P_VALS <- seq(0,1,1e-3)
RO_VALS <- seq(0, Phi, 0.0025)
CEX_MAIN <- 2.5
CEX_LAB <- 1.4

fixedvals <- function(x) { return (x[!is.infinite(x) & !is.na(x)]) }

calc.z0 <- function(rho, p)
{
  gz <- c( (1-p)*rho^2, #z^3
              -(rho^2 + (3-2*p)*rho), #z^2
                (2*rho + 2), #z
                  -1 ) #1
  roots <- Re(polyroot(rev(gz)))
  z0 <- roots[roots>0 & roots<1]
  stopifnot(length(z0) == 1)

  return(z0)
}

calc.z0.p1 <- function(rho)
{
  sqrrp1 <- sqrt(rho+1)
  d <- rho*sqrrp1
  
  return ( (sqrrp1-1)/d )
}

p0dot <-function(rho, p) 1 / (1+p*rho)
p1dot <-function(rho, p) 1 - p0dot(rho,p)
rho.hat <- function(rho, p) rho - p1dot(rho,p)
p.min <- function(rho) (rho-1)/(rho*(2-rho))

calc.p00 <- function(rho, p, z0=NULL)
{
  if (is.null(z0))
    z0 <- calc.z0(rho, p)
  
  rho.h <- rho.hat(rho,p)
  p00 <- (1-rho.h)*z0 / ( (1-z0)*(1-rho*z0) )
  stopifnot(p00 >= 0 & p00 <= 1)
  return (p00)
}

calc.p00.p1 <-function(rho)
{
  g <- 1+rho-rho^2
  
  sqrrp1 <- sqrt(rho+1)
  
  f <- (sqrrp1-1) / (rho^2-1+sqrrp1)
  
  return (f*g)
}

calc.p01 <- function(rho, p, z0=NULL)
{
  p00 <- calc.p00(rho, p, z0)
  
  p01 <- 1 - rho.hat(rho, p) - p00
  stopifnot(p01 >= 0 & p01 <= 1)
  return (p01)
}

calc.e<- function(rho, p, z0=NULL)
{
  if (p <= p.min(rho))
    return (Inf)
  
  if (p == 0)
    return (rho / (1-rho))
  
  if (is.null(z0))
    z0 <- calc.z0(rho, p)
  
  p0 <- p0dot(rho,p)
  p1 <- p1dot(rho,p)
  e0 <- calc.e0(rho,p,z0)
  e1 <- calc.e1(rho,p,z0)
  e <- p0*e0 + p1*e1
  
  if (is.valid.positive(e))
    return (e)
  
  return (NA)
}

calc.e0 <- function(rho, p, z0=NULL)
{
  if (p <= p.min(rho))
    return (Inf)
  
  if (p == 0)
    return (rho / (1-rho))
  
  if (is.null(z0))
    z0 <- calc.z0(rho, p)
  
  rho.h <- rho.hat(rho,p)
  p00 <- calc.p00(rho, p, z0)
  
  g1 <- (1-p)*rho^2 - (rho^2+(3-2*p)*rho) + 2*rho + 1
  gPrime1 <- 3*(1-p)*rho^2 - 2*(rho^2+(3-2*p)*rho) + (2*rho+2)
  G1 <- 1-rho.h
  GPrime1 <- G1 + (1-rho)*p00
  
  e0 <- 1/p0dot(rho,p) * (GPrime1*g1 - gPrime1*G1) / ((g1)^2)
  
  if (is.valid.positive(e0))
    return (e0)
  
  return (e0)
}

calc.e0.p1 <- function(rho)
{
  sqrrp1 <- sqrt(rho+1)
  
  return ( rho*(sqrrp1-1)/(sqrrp1-rho) )
}

is.valid.positive <- function(t) !is.null(t) & !is.na(t) & is.finite(t) & t >= 0

calc.e1 <- function(rho, p, z0=NULL)
{
  if (p <= p.min(rho))
    return (Inf)
  
  if (p == 0)
    return (NA)
  
  if (is.null(z0))
    z0 <- calc.z0(rho, p)
  
  rho.h <- rho.hat(rho,p)
  p01 <- calc.p01(rho, p, z0)
  
  g1 <- (1-p)*rho^2 - (rho^2+(3-2*p)*rho) + 2*rho + 1
  gPrime1 <- 3*(1-p)*rho^2 - 2*(rho^2+(3-2*p)*rho) + (2*rho+2)
  G1 <- p*rho*(1-rho.h)
  GPrime1 <- G1 + (1-(1-p)*rho)*p01
  
  e1 <- 1/p1dot(rho,p) * (GPrime1*g1 - gPrime1*G1) / ((g1)^2)
  
  if (is.valid.positive(e1))
    return (e1)
  
  return (NA)
}

calc.e1.rho1 <- function(p)
{
  p01 <- calc.p01(rho = 1, p = p)
  return (((1+p)*p01 + 1)/p )
}

calc.p.eq.rho1 <- function(gamma) min((sqrt(1+4*gamma)-1)/2, 1)

plot.queue.length <- function(rho, reso = 0.005, xlim = c(0, 1), ylim = NULL)
{
  f1 <- function(p) calc.e(rho, p)
  f2 <- function(p) calc.e0(rho, p)
  f3 <- function(p) calc.e1(rho, p)
  f4 <- function(p) mm1.approximation(rho, p)
  
  p <- seq(0, 1, reso)
  Lp1 <- sapply(p, f1)
  Lp2 <- sapply(p, f2)
  Lp3 <- sapply(p, f3)
  Lp4 <- sapply(p, f4)
  
  if(is.null(ylim))
    ylim <- c(min(fixedvals(c(Lp1, Lp2))), max(fixedvals(c(Lp1, Lp2))))
  
  colors <- c("red","blue", "black")
 
  plot(p, Lp1, type='l', ylim = ylim, xlim = xlim, ylab = "")
  par(new=T)
  plot(p, Lp2, type='l', col = "blue", ylim = ylim, xlim = xlim, ylab = "")
  par(new=T)
  plot(p, Lp3, type='l', col = "red", ylim = ylim, xlim = xlim, ylab = "")
  par(new=T)
  plot(p, Lp4, type='l', col = "gray", ylim = ylim, xlim = xlim, ylab = "", lty = "dotted")
  abline(v = ((rho - 1) / (rho * (2-rho))))
 
 legend("topright", c("E[L | Y=1]", "E[L | Y=0]", "E[L]"), col=colors, text.col=colors, y.intersp = 0.2, bty='n')
 
  
  names(Lp2) <- p
  return (Lp2)
}

calc.p.eq <- function(rho, gamma)
{
  if (rho <= 1 / (1 + gamma))
    return (0)
  
  if (gamma*calc.e0(rho, 1) >= 1+rho)
    return (1)
  
  f <- function(p)
  {
    e0 <- calc.e0(rho, p)
    return (abs(gamma*e0/(1+p*rho) - 1))
  }
  
  return (optimize(f, c(0,1))$minimum)
}

calc.p.eq.vec <- function(gamma)
{
  f <- function(rho) {return (calc.p.eq(rho, gamma))}
  return(sapply(RO_VALS, f))
}

plot.costs.chart <- function(rho, gamma, reso = 0.005)
{
  f <- function(p) { return (calc.e(rho, p)) }
  f0 <- function(p) { return (calc.e0(rho, p)) }
  
  p <- seq(0, 1, reso)
  Lp <- sapply(p, f)
  Lp0 <- sapply(p, f0)
  
  Vs_p <- (1 + gamma*Lp0*p*rho/(1+p*rho))
  Vn_p <- gamma*Lp
  
  ylimLo <- min(fixedvals(c(Vs_p, Vn_p)))
  ylimUp <- max(fixedvals(c(Vs_p, Vn_p)))
  
  if (rho >= 1)
    ylimUp <- min(ylimUp, 15*ylimLo)
  

  d <- (ylimUp-ylimLo)/20
  ylimLo <- ylimLo - d
  
  if (tail(fixedvals(Vs_p),1) < tail(fixedvals(Vn_p),1))
    d <- -d
  
  plot(p, Vs_p, col="gray70", type='l', xlab="", ylab="", xlim=c(0,1), ylim=c(ylimLo, ylimUp), lwd = 2)
  
  text(x = 0.95, y = tail(fixedvals(Vs_p),1)+d, labels = bquote(C[S]), cex = 1.7, col = "gray70")
  par(new=T)
  plot(p, Vn_p, col="gray10", type='l', xlab="", ylab="", xlim=c(0,1), ylim=c(ylimLo, ylimUp), lwd = 2)
  text(x = 0.95, y = tail(fixedvals(Vn_p),1)-d, labels = bquote(C[N]), cex = 1.7, col = "gray10")
  title(main = bquote(C[N] ~ " & " ~ C[S] ~ " v.s." ~ p), cex.main = CEX_MAIN,
        ylab = "Cost", cex.lab = CEX_LAB, xlab = bquote(p))
}

require('plotrix')

accumulated.p.ro.chart <- function()
{
  gamma_vals <- sort(c(1, 2, 5, 10, .1, .2, .5 ), decreasing = F)
  ngraphs <- length(gamma_vals)
  colors <- c("gray0", "gray45", "gray70")
  i<-1
  
  for (gamma in gamma_vals)
  {
    p_eqs <- calc.p.eq.vec(gamma)
    color <- colors[(i%%length(colors))+1]
    
    plot(RO_VALS, p_eqs, xlab = "", ylab = "", xlim = c(0, Phi), ylim = c(0,1), type = 'l', lwd = 2, pch = 20, col = color)
    
    txtanchor <- which(p_eqs >= 0.1*ngraphs-0.075*(ngraphs-i))[1]
    tx <- RO_VALS[txtanchor]
    ty <- p_eqs[txtanchor]

    text(x = tx, y = ty+0.05, labels = bquote(gamma  == ~ .(gamma)), cex = 1.2, col = color)
    i <- i+1
    par(new = T)
  }

  title(main = bquote(p[e] ~ v.s. ~ rho), cex.main = CEX_MAIN, 
        xlab = bquote(rho), ylab = bquote(p[e]), cex.lab = CEX_LAB)
  par(new = F)
  
  p1s <- sapply(X = gamma_vals, calc.first.pe1.rho)
  points(p1s, rep(1,length(gamma_vals)), col='red')
}

average.cost <- function(rho, gamma, p)
{
  if(rho.hat(rho, p) >= 1)
    return (Inf)
  
  e <- calc.e(rho, p)
  e0 <- calc.e0(rho, p)
  
  return (p/gamma + e - p/(1+p*rho)*e0)
}

plot.tmp <- function (rho)
{
  f <- function(p)
  {
    if(rho.hat(rho, p) >= 1)
      return (Inf)
    
    e <- calc.e(rho, p)
    e0 <- calc.e0(rho, p)
    
    return (e - p/(1+p*rho)*e0)
  }
  
  plot(P_VALS, sapply(P_VALS, f), type='l')
}

approximated.average.cost <- function(rho, gamma, p)
{
  rho.h <- rho.hat(rho, p)
  if(rho.h >= 1)
    return (Inf)
#   em <- mm1.approximation(rho, p)
#   e <- mm1.queue.length((1-p)*rho)
#   p1 <- p1dot(rho, p)
  pr0 <- p0dot(rho, p)
#   e0 <- calc.e0(rho, p)
#   e1 <- mm1.queue.length(rho)
#   pe <- calc.p.eq(rho, gamma)
  return (p/gamma + (1-p*pr0)*rho.h/(1-rho.h))
#   return (p  + gamma*rho/(1+rho)*e1)
}

calc.social.opt <- function(rho, gamma)
{
  if (rho == 0)
    return (0)
  
  if (gamma == 0)
    return (0)
  
  f <- function(p) average.cost(rho, gamma, p)
 
  return (optimize(f, c(0, 1))$minimum)
}

calc.approximated.social.opt <- function(rho, gamma)
{
  if (rho == 0)
    return (0)
  
  if (gamma == 0)
    return (0)
  
  f <- function(p) approximated.average.cost(rho, gamma, p)
  
  return (optimize(f, c(0, 1))$minimum)
}

plot.average.cost <- function(gamma, rho, xlim=NULL, ylim=NULL)
{
  minarg <- calc.social.opt(rho, gamma)
  minval <- average.cost(rho, gamma, minarg)
  f <- function(p) average.cost(rho, gamma, p)
  
  if (is.null(xlim))
    xlim<-c(0,1)
  
  if (is.null(ylim))
    ylim <- c(minval,f(0))
  
  plot(P_VALS, sapply(P_VALS, f), ylim = ylim, xlim=xlim, type='p')
  points(minarg, minval, col='red', lwd =3)  
  text(x = minarg, y = minval, labels = paste("(", round(minarg,3), "; ", round(minval,3), ")", sep = ""), 
       cex = 1, col = "red", pos = 3 )
  grid()
  
  par(new=T)
  plot(P_VALS, sapply(P_VALS, function(p) return(calc.e(rho,p))), ylim = ylim, xlim=xlim, type='p', xlab = "", ylab = "")
}



plot.approximated.average.cost <- function(gamma, rho)
{
  minarg <- calc.approximated.social.opt(rho, gamma)
  minval <- approximated.average.cost(rho, gamma, minarg)
  f <- function(p) approximated.average.cost(rho, gamma, p)
  plot(P_VALS, sapply(P_VALS, f), xlim = c(0,1), type='l')
  points(minarg, minval, col='red', lwd =3)  
  text(x = minarg, y = minval, labels = paste("(", round(minarg,3), "; ", round(minval,3), ")", sep = ""), 
       cex = 1, col = "red", pos = 3 )
  grid()
}

calc.social.opt.vec <- function(gamma)
{
  f <- function(rho) {return (calc.social.opt(rho, gamma))}
  
  return (sapply(RO_VALS, f))
}

calc.approximated.social.opt.vec <- function(gamma)
{
  f <- function(rho) {return (calc.approximated.social.opt(rho, gamma))}
  
  return (sapply(RO_VALS, f))
}


is.p1.eq <- function(rho, gamma)
{
  theta <- sqrt(1+rho)
  
  x <- theta^4 + (gamma-1)*theta^3 - (gamma+1)*theta^2 - gamma*theta + gamma
  
  return(x)
}

calc.first.pe1.rho <- function(gamma)
{
  P <- c(1, gamma-1, -(gamma+1), -gamma, gamma)
  roots <- Re(polyroot(rev(P)))
  theta <- roots[roots>1 & roots<Phi]
  
  return(theta^2-1)
}

calc.zero.priority.waiting1 <- function(rho, p)
{
  e1 <- calc.e1(rho, p)
  rho_hat <- rho.hat(rho, p)
  
  return ((e1+1)/rho_hat - 1)
}

calc.zero.priority.waiting <- function(rho, p)
{
  e <- calc.e(rho, p)
  rho_hat <- rho.hat(rho, p)
  
  return ((e+1)/rho_hat - 1)
}

calc.tmp <- function(gamma, rho, p=1)
{
  w <- calc.zero.priority.waiting(rho, p)
  w1 <- calc.zero.priority.waiting1(rho, p)
  p0 <- p0dot(rho, p)
  
  v <- 1+gamma*((1-p0^2)*w1 - w)
  return (v)
}

plot.social.vs.individual.opt <- function(gamma)
{
  peqs <- calc.p.eq.vec(gamma)
  psos <- calc.social.opt.vec(gamma)
  psosapp <- calc.approximated.social.opt.vec(gamma)
  
  txpeqs <- RO_VALS[which(peqs >= 0.5)[1]]
  typeqs <- peqs[which(peqs >= 0.5)[1]]

  txpsos <- RO_VALS[which(psos >= 0.5)[1]]
  typsos <- psos[which(psos >= 0.5)[1]]
  
  plot(RO_VALS, peqs, ylim = c(0,1), type='l', lwd = 2, col = "gray70", ylab ="", xlab="")
  text(x = txpeqs, y = typeqs-.09, labels = bquote(p[e]), cex = 1.7, col = "gray70")
  
  par(new=T)
  
  plot(RO_VALS, psos, ylim = c(0,1), type='l', lwd = 2, col = "gray10", ylab ="", xlab="")
  text(x = txpsos, y = typsos+.05, labels = "p*", cex = 1.7, col = "gray10")

  par(new=T)  
  plot(RO_VALS, psosapp, ylim = c(0,1), type='l', lwd = 2, col = "red", lty = "dotted", ylab ="", xlab="")
  
  points(y=0, x=1-sqrt(gamma/(1+gamma)), col="red")

  title(main = bquote(p[e] ~ "&" ~ "p*" ~ " v.s." ~ rho), cex.main = CEX_MAIN,
        ylab = "Probability", cex.lab = CEX_LAB, xlab = bquote(rho))
}

calc.first.pe1.gamma <- function(rho)
{
  t <- sqrt(1+rho)
  
  res <- (t^2 + t^3 - t^4) / (1 - t - t^2 + t^3)
  return (res)
}

plot.rev.social.vs.individual.opt <- function(rho)
{
  GBOUND <- 1.2 * calc.first.pe1.gamma(rho)
  GAMMA_VALS <- seq(0, GBOUND, 0.005)
  
  if(rho > 1)
    GAMMA_VALS <- GAMMA_VALS[-1]
  
  fpeq <- function(gamma) calc.p.eq(rho, gamma)
  fpso <- function(gamma) calc.social.opt(rho, gamma)
  fpsoapp <- function(gamma) calc.approximated.social.opt(rho, gamma)
  
  peqs <- sapply(GAMMA_VALS, fpeq)
  psos <- sapply(GAMMA_VALS, fpso)
  psosapp <- sapply(GAMMA_VALS, fpsoapp)
  
  txpeqs <- GAMMA_VALS[which(peqs >= 0.5)[1]]
  typeqs <- peqs[which(peqs >= 0.5)[1]]
  
  txpsos <- GAMMA_VALS[which(psos >= 0.5)[1]]
  typsos <- psos[which(psos >= 0.5)[1]]
  
  plot(GAMMA_VALS, peqs, ylim = c(0,1), type='l', lwd = 2, col = "gray70", ylab ="", xlab="")
  text(x = txpeqs, y = typeqs-.09, labels = bquote(p[e]), cex = 1.7, col = "gray70")
  
  par(new=T)
  
  plot(GAMMA_VALS, psos, ylim = c(0,1), type='l', lwd = 2, col = "gray10", ylab ="", xlab="")
  text(x = txpsos, y = typsos+.05, labels = "p*", cex = 1.7, col = "gray10")
  
  title(main = bquote(p[e] ~ "&" ~ "p*" ~ " v.s." ~ gamma), cex.main = CEX_MAIN,
      ylab = "Probability", cex.lab = CEX_LAB, xlab = bquote(gamma))

  h <- 1e-6
  
  f <- function(p) calc.e(rho, p) - p*p0dot(rho, p)*calc.e0(rho, p)
  
  gamma1appval <- h / (f(1-h) - f(1))
  
  abline(v=gamma1appval, lty="dotted")
}

generate.costs.chart <- function (gamma, rho)
{
  filename <- chartr('.','_', paste("cost_vs_p_", gamma, "_", rho, sep=""))
  jpeg(paste("plots/", filename, ".png" , sep=""), width = 1000, height = 1200)
  plot.costs.chart(rho, gamma)
  dev.off()
}

generate.pe.vs.pstar.chart <- function (gamma)
{
  filename <- chartr('.','_', paste("pe_vs_pstar_", gamma, sep=""))
  jpeg(paste("plots/", filename, ".png" , sep=""), width = 1000, height = 1200)
  plot.social.vs.individual.opt(gamma)
  dev.off()
}

generate.pe.vs.pstar.rev.chart <- function (rho)
{
  filename <- chartr('.','_', paste("pe_vs_pstar_rev_", rho, sep=""))
  jpeg(paste("plots/", filename, ".png" , sep=""), width = 1000, height = 1200)
  plot.rev.social.vs.individual.opt(rho)
  dev.off()
}

calc.p1e1 <- function(rho, p) calc.e1(rho, p)*p1dot(rho, p)

mm1.queue.length <- function(rho) rho/(1-rho)
mm1.approximation <- function(rho, p) mm1.queue.length(rho.hat(rho, p))

calc.p1e1.vec <- function(rho)
{
  return (sapply(P_VALS, function(p) calc.p1e1(rho, p)))
}