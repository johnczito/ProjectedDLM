# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

set.seed(8675309)

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

datasource = "Black Mountain"

if(datasource == "Black Mountain"){
  raw_data <- read.csv("datasets/black_mountain_wind_direction.csv", header = FALSE)
  my_radians_0_2pi <- degrees2radians(raw_data$V1)
} else if (datasource == "O'Hare"){
  raw_data <- read.csv("datasets/ohare.csv", header = TRUE) |>
    fill(HourlyWindDirection)
  my_radians_0_2pi <- degrees2radians(raw_data$HourlyWindDirection)
}

U <- radians2unitcircle(my_radians_0_2pi)

# ------------------------------------------------------------------------------
# basic settings
# ------------------------------------------------------------------------------

n = ncol(U)
T = nrow(U)
p = n
FF = getIDFF(T, n)

# ------------------------------------------------------------------------------
# PDLM settings
# ------------------------------------------------------------------------------

FF = getIDFF(T, n)

pdlm_ndraw = 5000
pdlm_burn  = 5000
pdlm_thin  = 10

pdlm_draws = gibbs_pdlm(U, FF, 
                        ndraw = pdlm_ndraw, 
                        burn = pdlm_burn, 
                        thin = pdlm_thin)

par(mfrow = c(n, n))

for(i in 1:n){
  for(j in 1:n){
    plot(1:pdlm_ndraw, pdlm_draws$G[i, j, ], type = "l")
  }
}

par(mfrow = c(n, n))

for(i in 1:n){
  for(j in 1:n){
    plot(1:pdlm_ndraw, pdlm_draws$W[i, j, ], type = "l")
  }
}

par(mfrow = c(n, n))

for(i in 1:n){
  for(j in 1:n){
    plot(1:pdlm_ndraw, pdlm_draws$Sigma[i, j, ], type = "l")
  }
}
