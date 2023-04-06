myqqplot <- function(x, y, p = seq(0.05, 0.995, length.out = 250),
                     xlab = "", ylab = "", main = "", pch = 19, cex = 0.25){
  qx = quantile(x, probs = p)
  qy = quantile(y, probs = p)
  plot(qx, qy, xlab = xlab, ylab = ylab, pch = pch, cex = cex, main = main)
  abline(a = 0, b = 1, col = "red", lty = 2)
}

# WIPWIPWIP

plot_circular_band <- function(xpts, upper, lower, color){
  TT = length(upper)
  new_upper = upper
  upper2 = rep(0, TT)
  lower2 = rep(0, TT)
  for (t in 1:TT){
    if (lower[t] > upper[t]){
      new_upper[t] = 2 * pi
      upper2[t] = upper[t]
      lower2[t] = 0
    }
  }
  polygon(c(xpts, rev(xpts)), c(new_upper, rev(lower)),
          col = color, lty = 0)
  polygon(c(xpts, rev(xpts)), c(upper2, rev(lower2)),
          col = color, lty = 0)
}