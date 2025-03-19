# ==============================================================================
# load packages
# ==============================================================================

library(Ternary)

# ==============================================================================
# load data
# ==============================================================================

Energy <- read.csv("datasets/annual_generation_state_1990_2022.csv", header = TRUE)
energy <- as.data.frame(Energy)
index <- which(energy$`ENERGY.SOURCE`=="Total")
energy <- energy[-index,]

# ==============================================================================
# data pre-processing
# ==============================================================================

years <- sort(unique(energy$YEAR))

W <- sapply(years, function(yr){
  
  index <- which(energy==yr)
  temp <- energy[index,]
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Coal", "Petroleum"))
  coal <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Natural Gas"))
  gas <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Natural Gas", "Coal", "Petroleum"))
  others <- sum(temp$`GENERATION..Megawatthours`[-index], na.rm = TRUE)
  
  re <- c(coal, gas, others)
  return(re/sum(re))
})

T = ncol(W)
n = nrow(W)

labs = c("coal/petroleum", "natural gas", "nuclear/renewables")
cls = c("red", "darkgreen", "blue")

# ==============================================================================
# plot the shares
# ==============================================================================

par(mfrow = c(1, 2), mar = c(2, 4, 2, 2))

plot(years, W[1, ], type = "l", ylim = c(0, 1), ylab = "Share", xlab = "", 
     col = cls[1], lwd = 2)

for(i in 2:n){
  lines(years, W[i, ], col = cls[i], lwd = 2)
}

legend("topright", labs, col = cls, lwd = 1, bty = "n", lty = 1)

par(mar = c(0, 0, 0, 0))

TernaryPlot(alab = paste(labs[1], "\u2192"), 
            blab = paste(labs[2], "\u2192"), 
            clab = paste("\u2190", labs[3]),
            lab.col = cls,
            point = "up", 
            lab.cex = 0.8, grid.minor.lines = 0,
            grid.lty = "solid", col = rgb(0.9, 0.9, 0.9), grid.col = "white", 
            axis.col = rgb(0.6, 0.6, 0.6), ticks.col = rgb(0.6, 0.6, 0.6),
            axis.rotate = FALSE,
            padding = 0.08)
cols <- TernaryPointValues(rgb)
ColourTernary(cols, spectrum = NULL)
AddToTernary(graphics::points, 100*t(W), pch = 19, col = "white", cex = 0.5)
TernaryArrows(100*t(W)[T, ] + c(0, 15, -15), 100*t(W)[T, ] - c(0, -5, 5), length = 0.1, col = "white")
TernaryArrows(100*t(W)[1, ] + c(0, 15, -15), 100*t(W)[1, ] - c(0, -5, 5), length = 0.1, col = "white")
AddToTernary(text, 100*t(W)[1, ] + c(0, 22, -22), years[1], cex = 0.75, font = 2, col = "white")
AddToTernary(text, 100*t(W)[T, ] + c(0, 22, -22), years[T], cex = 0.75, font = 2, col = "white")


