# prob.as.a.function.of.p00.and.p01 <- function(p, ro, n)
# {
#   sumlist <- function(l) { return(Reduce("+",l)) }
#   lp_n0 <- list(c(1,0))
#   lp_n1 <- list(c(0,1))
#   
#   if (n > 0)
#   {
#     for(i in 2:(n+1))
#     {
#       s <- sumlist(head(lp_n1, i-1)) - ro*p*sumlist(head(lp_n0, i-1))
#       lp_n1[[i]] <- ro*lp_n1[[i-1]] + s
#       lp_n0[[i]] <- (1-p)*ro*lp_n0[[i-1]] - s
#     }
#   }
#   
#   cat("P_",n,"0 = ", sep="")
#   show(lp_n0[[n+1]])
#   cat("P_",n,"1 = ", sep="")
#   show(lp_n1[[n+1]])
#   
#   return (lp_n0[[n+1]]+lp_n1[[n+1]])
# }
# 
# p11.and.p10 <- function(p, ro)
# {
#   return ( ro*c(1-p, 1) )
# }
# 
# p21.and.p20 <- function(p, ro)
# {
#   t1 <- ro - 2*p*ro
#   t2 <- p + ro
#   return ( ro*c(t1, t2) )
# }
# 
# p31.and.p30 <- function(p, ro)
# {
#   t1 <- ro^2 - 3*p*ro^2 - 2*p^2*ro
#   t2 <- 3*p*ro + 2*p + ro^2
#   return ( ro*c(t1, t2) )
# }
# 
# p41.and.p40 <- function(p, ro)
# {
#   t1 <- ro^3 - 4*p*ro^3 - 7*p^2*ro^2 - 4*p^2*ro
#   t2 <- 7*p*ro + 2*p^2*ro + 6*p*ro^2 + 4*p + ro^3
#   return ( ro*c(t1, t2) )
# }
# 
# p51.and.p50 <- function(p, ro)
# {
#   t1 <- ro^4 - 16*p^2*ro^2 - 5*p*ro^4 - 16*p^2*ro^3 - 4*p^3*ro^2 - 8*p^2*ro
#   t2 <- 16*p*ro^2 + 16*p*ro + ro^4 + 10*p*ro^3 + 9*p^2*ro^2 + 8*p^2*ro + 8*p
#   return ( ro*c(t1, t2) )
# }