# ------------------------------------------------------------------------------
# packages
# ------------------------------------------------------------------------------

source("_packages.R")

# https://plotly.com/r/3d-subplots/

# ------------------------------------------------------------------------------
# function that evaluates the pdf of the projected normal in the n = 2 case.
# the formula comes from Wang, Gelfand (2013 SM)
# ------------------------------------------------------------------------------

dpn <- function(o, m1, m2, v1, v2, r){
  s1 <- sqrt(v1)
  s2 <- sqrt(v2)
  m  <- c(m1, m2)
  V  <- cbind(c(v1, r * s1 * s2), c(r * s1 * s2, v2))
  a  <- 1 / (s1 * s2 * sqrt(1 - r^2))
  C  <- (a^2) * (v2 * cos(o)^2 - r * s1 * s2 * sin(2*o) + v1 * sin(o)^2)
  D  <- (a^2) * (m1 * s2 * (s2 * cos(o) - r * s1 * sin(o)) + m2 * s1 * (s1 * sin(o) - r * s2 * cos(o))) / sqrt(C)
  k  <- dmvnorm(m, sigma = V) + a * D * pnorm(D) * dnorm(a * (m1 * sin(o) - m2 * cos(o)) / sqrt(C))
  f <- k / C
  return(f)
}

# ------------------------------------------------------------------------------
# boiler plate
# ------------------------------------------------------------------------------

n <- 500
o <- seq(0, 2*pi, length.out = n)
x <- cos(o)
y <- sin(o)
z1 <- numeric(n)

# ------------------------------------------------------------------------------
# Panel 1: unimodal, asymmetric example
# ------------------------------------------------------------------------------

m1 = 2
m2 = 0
v1 = 1
v2 = 1
r  = 0.75

z2 <- dpn(o, m1, m2, v1, v2, r)

data <- data.frame(x, y, z1, z2)

fig1 <- plot_ly(data, x = ~x, y = ~y, z = ~z2, type = 'scatter3d', mode = 'lines',
                line = list(width = 6, color = "red"))
fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F),
                                     yaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F),
                                     zaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F)))

fig1 <- fig1 %>% add_trace(x = ~x, y = ~y, z = ~z1, data = data,
                           line = list(width = 6, color = "black"))


fig1

# ------------------------------------------------------------------------------
# Panel 2: antipodal, but unequal
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Panel 3: bimodal, not antipodal
# ------------------------------------------------------------------------------

m1 = -0.25
m2 = -0.25
v1 = 1
v2 = 1
r  = -0.75

z2 <- dpn(o, m1, m2, v1, v2, r)

data <- data.frame(x, y, z1, z2)

fig1 <- plot_ly(data, x = ~x, y = ~y, z = ~z2, type = 'scatter3d', mode = 'lines',
                line = list(width = 6, color = "red"))
fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F),
                                     yaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F),
                                     zaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F)))

fig1 <- fig1 %>% add_trace(x = ~x, y = ~y, z = ~z1, data = data,
                           line = list(width = 6, color = "black"))


fig1

# ------------------------------------------------------------------------------
# Panel 4: aymmetric, nearly bimodal
# ------------------------------------------------------------------------------

m1 = -0.1
m2 = -0.1
v1 = 0.015
v2 = 0.01
r  = -0.5

z2 <- dpn(o, m1, m2, v1, v2, r)

data <- data.frame(x, y, z1, z2)

fig1 <- plot_ly(data, x = ~x, y = ~y, z = ~z2, type = 'scatter3d', mode = 'lines',
                line = list(width = 6, color = "red"))
fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F),
                                     yaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F),
                                     zaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F)))

fig1 <- fig1 %>% add_trace(x = ~x, y = ~y, z = ~z1, data = data,
                           line = list(width = 6, color = "black"))


fig1

# ------------------------------------------------------------------------------
# Panel 5: antipodal, diffuse (potato chip)
# ------------------------------------------------------------------------------

m1 = 0
m2 = 0
v1 = 1
v2 = 1
r  = 0.5

z2 <- dpn(o, m1, m2, v1, v2, r)

data <- data.frame(x, y, z1, z2)

fig1 <- plot_ly(data, x = ~x, y = ~y, z = ~z2, type = 'scatter3d', mode = 'lines',
                line = list(width = 6, color = "red"))
fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F),
                                     yaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F),
                                     zaxis = list(title = '', showgrid = F, zeroline = F, showticklabels = F)))

fig1 <- fig1 %>% add_trace(x = ~x, y = ~y, z = ~z1, data = data,
                           line = list(width = 6, color = "black"))


fig1
