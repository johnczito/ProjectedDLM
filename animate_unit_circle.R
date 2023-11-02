# https://www.datanovia.com/en/blog/gganimate-how-to-create-plots-with-beautiful-animation-in-r/
# https://towardsdatascience.com/how-to-create-animated-plots-in-r-adf53a775961

source("helpers/_helpers.R")

library(animation)

black_mountain_degrees <- read.csv("datasets/black_mountain_wind_direction.csv", header = FALSE)
black_mountain_radians <- degrees2radians(black_mountain_degrees$V1)
black_mountain_unit_circle <- radians2unitcircle(black_mountain_radians)
U <- black_mountain_unit_circle

resolution = 10

angles = c(seq(pi / 4, 0, length.out = resolution),
           seq(2*pi, 0, length.out = 8 * resolution),
           seq(2*pi, 7*pi/4, length.out = resolution))
angles = sample(c(runif(90, 7*pi/4, 2*pi), runif(10, 0, pi/4)),
                size = 100,
                replace = TRUE)
#angles = black_mountain_radians
TT = length(angles)

saveGIF(
  {

    for(t in 1:TT){

      par(mfrow = c(1, 2))

      # plot unit circle

      unit_circle_x_vals = seq(-1, 1, length.out = 100)
      unit_circle_y_vals = sqrt(1 - unit_circle_x_vals ^ 2)

      plot(c(unit_circle_x_vals, rev(unit_circle_x_vals)),
           c(unit_circle_y_vals, -unit_circle_y_vals),
           type = "l",
           main = expression(paste(u[t], " = [cos", theta[t], " sin", theta[t], "]'")),
           xaxt = "n",
           xlab = "",
           yaxt = "n",
           ylab = "",
           bty = "n",
           cex.main = 2)
      grid()
      abline(v = 0, h = 0)
      mtext(expression(3*pi/2), side = 1, line = 0)
      mtext(expression(pi), side = 2, line = 0)
      mtext(expression(pi/2), side = 3, line = 0)
      mtext(0, side = 4, line = 0)

      angle = angles[t]
      points(cos(angle), sin(angle), pch = 19, col = "red")
      segments(0, 0, cos(angle), sin(angle), col = "red", lty = 2)

      # plot time series of angles

      plot(1:t, angles[1:t], col = "red", type = "l",
           main = expression(paste(theta[t], " = atan2(", u[t], ") mod 2", pi)),
           ylab = "angle",
           xlab = "time",
           xlim = c(1, TT),
           ylim = c(0, 2*pi),
           yaxt = "n",
           cex.main = 2)
      axis(2, at = c(0, pi/2, pi, 3*pi/2, 2*pi),
           labels = c("0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))

    }
  },
  movie.name = "test.gif",
  interval = 0.1,
  ani.width = 1000 * 0.75,
  ani.height = 625 * 0.75,
  outdir = getwd()
)

# fine tune my gif
# generate preliminary interval results
# write section of paper about old methods
# write section of paper about forecasting
