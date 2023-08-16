# ------------------------------------------------------------------------------
# packages and helpers
# ------------------------------------------------------------------------------

library(circular)

change_0_2pi_to_minpi_pi <- function(x){
  n = length(x)
  y = x
  y[x > pi] = x[x > pi] - 2*pi
  return(y)
}

# ------------------------------------------------------------------------------
# fake data
# ------------------------------------------------------------------------------

set.seed(8675309)

c = circular( c(rvonmises(4500, 1.75 * pi, 2), rvonmises(500, 2.5, 6)), modulo = "2pi")

# ------------------------------------------------------------------------------
# quantiles
# ------------------------------------------------------------------------------

L = 0.1
U = 0.9
my_breaks = 75

my_qL = quantile.circular(c, probs = L)
my_median = median.circular(c)
my_qU = quantile.circular(c, probs = U)

# ------------------------------------------------------------------------------
# Panel 1: original sample and its circular quantiles
# ------------------------------------------------------------------------------

hist(as.vector(c), breaks = my_breaks, freq = FALSE, xlab = "x", main = "(A) Original sample", xaxt = "n", xlim = c(-pi, 2*pi))
abline(v = c(my_qL, my_median, my_qU), col = c("red", "black", "blue"), lwd = 3)
axis(1, at = c(-pi, -pi/2, 0, pi/2, pi, 3*pi/2, 2*pi),
     labels = c(expression(-pi), expression(-pi/2), "0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))

# ------------------------------------------------------------------------------
# Panel 2: rearrange sample so branch cut occurs at circular median
# ------------------------------------------------------------------------------

x = as.vector(c)
circularmedian = as.vector(my_median)
tx <- (x-circularmedian)%%(2*pi)

hist(tx, breaks = my_breaks, freq = FALSE, xlab = "y", main = "(B) Branch cut at median",  xaxt = "n", xlim = c(-pi, 2*pi))
axis(1, at = c(-pi, -pi/2, 0, pi/2, pi, 3*pi/2, 2*pi),
     labels = c(expression(-pi), expression(-pi/2), "0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))


# ------------------------------------------------------------------------------
# Panel 3: "linearize"
# ------------------------------------------------------------------------------

tx <- change_0_2pi_to_minpi_pi(tx)

hist(tx, breaks = my_breaks, freq = FALSE, xlab = "y", main = '(C) "Linearized" sample',  xaxt = "n", xlim = c(-pi, 2*pi))
abline(v = quantile(tx, probs = c(L, 0.5, U)), col = c("red", "black", "blue"), lwd = 3)
axis(1, at = c(-pi, -pi/2, 0, pi/2, pi, 3*pi/2, 2*pi),
     labels = c(expression(-pi), expression(-pi/2), "0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))

# ------------------------------------------------------------------------------
# Compare
# ------------------------------------------------------------------------------

linearized_quantiles <- quantile(tx, probs = c(L, U))
my_quantile = (linearized_quantiles + circularmedian)%%(2*pi)

# THESE SHOULD MATCH

my_quantile
c(my_qL, my_qU)
