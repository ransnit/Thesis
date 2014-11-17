plot.best.response <- function(p)
{
  id <- function(x){return(x)}
  plot.function(id, from = 0, to = 1, lty = "dashed")
  segments(x0 = 0, y0 = 0, x1 = p, y1 = 0)
  segments(x0 = p, y0 = 0, x1 = p, y1 = 1)
  segments(x0 = p, y0 = 1, x1 = 1, y1 = 1)
}