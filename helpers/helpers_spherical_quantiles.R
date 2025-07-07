perpendicular_circle_intersections <- function(m, c) {
  if (abs(sqrt(sum(m^2)) - 1) > 1e-8) stop("m must be a unit vector")
  if (length(m) != 2) stop("m must be a 2D vector")
  if (abs(c) >= 1) stop("c must be between -1 and 1 (exclusive)")
  
  p <- c * m
  v <- c(-m[2], m[1])
  
  a <- sum(v^2)
  b <- 2 * sum(p * v)
  d <- sum(p^2) - 1
  
  disc <- b^2 - 4 * a * d
  if (disc < 0) return(NULL)
  
  t1 <- (-b + sqrt(disc)) / (2 * a)
  t2 <- (-b - sqrt(disc)) / (2 * a)
  
  pt1 <- p + t1 * v
  pt2 <- p + t2 * v
  
  # Use cross product to determine which is counterclockwise from m
  cross1 <- m[1] * pt1[2] - m[2] * pt1[1]
  cross2 <- m[1] * pt2[2] - m[2] * pt2[1]
  
  # Cross product > 0 means pt is counterclockwise from m
  if (cross1 > cross2) {
    rbind(pt1, pt2)  # pt1 is CCW from m
  } else {
    rbind(pt2, pt1)  # pt2 is CCW from m
  }
}

plot_circular_prediction_set <- function(draws, obs, alpha, main, ymax = NULL){
  angle_draws = unitcircle2radians(draws)
  med = mediandir(draws)
  c = quantile(c(draws %*% med), alpha)
  C <- perpendicular_circle_intersections(med, c)
  int_bounds <- unitcircle2radians(C)
  int_U <- int_bounds[1]
  int_L <- int_bounds[2]
  if(is.null(ymax)){
    hist(angle_draws, breaks = "Scott", freq = FALSE, main = main,
         border = "white")
    ymax <- par("usr")[4]
  }else{
    hist(angle_draws, breaks = "Scott", freq = FALSE, main = main,
         border = "white", ylim = c(0, ymax))
  }
  
  abline(v = my_radians_0_2pi[t], col = "blue", lwd = 2)
  ymin <- par("usr")[3]
  
  if(int_L < int_U){
    polygon(
      x = c(int_L, int_U, int_U, int_L),
      y = c(ymin, ymin, ymax, ymax),
      col = rgb(1, 0.5, 0, 0.25),
      border = NA
    )
    th <- seq(int_L, int_U, length.out = 200)
  }else{
    polygon(
      x = c(int_L, 2*pi, 2*pi, int_L),
      y = c(ymin, ymin, ymax, ymax),
      col = rgb(1, 0.5, 0, 0.25),
      border = NA
    )
    polygon(
      x = c(0, int_U, int_U, 0),
      y = c(ymin, ymin, ymax, ymax),
      col = rgb(1, 0.5, 0, 0.25),
      border = NA
    )
    th <- c(seq(int_L, 2*pi, length.out = 100), seq(0, int_U, length.out = 100))
  }
  
  arc <- radians2unitcircle(th)
  
  plot(draws, pch = 19, cex = 0.25, bty = "n")
  
  lines(arc, lwd = 3, col = rgb(1, 0.5, 0))
  arrows(0, 0, med[1], med[2], col = "red", lwd = 2)
  points(t(obs), col = "blue", pch = 19, lwd = 3)
}
