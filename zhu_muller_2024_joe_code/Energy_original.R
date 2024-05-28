library(readxl)
library(rgl)
library(maditr)
library(sphereplot)
library(manifold)

Energy <- read_excel("Energy.xls")
energy <- as.data.frame(Energy)
index <- which(energy$`ENERGY SOURCE`=="Total")
energy <- energy[-index,]

W <- sapply(c(2005:2019), function(yr){
  
  index <- which(energy==yr)
  temp <- energy[index,]
  
  index <- which(temp$`ENERGY SOURCE` %in% c("Coal", "Petroleum"))
  coal <- sum(temp$`GENERATION (Megawatthours)`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY SOURCE` %in% c("Natural Gas"))
  gas <- sum(temp$`GENERATION (Megawatthours)`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY SOURCE` %in% c("Natural Gas", "Coal", "Petroleum"))
  others <- sum(temp$`GENERATION (Megawatthours)`[-index], na.rm = TRUE)
  
  re <- c(coal, gas, others)
  return(re/sum(re))
})

truth <- sqrt(W[,dim(W)[2]])
Z <- sqrt(W[,c(1:(dim(W)[2]-1))])

mfd <- structure(1, class='Sphere')
mu <- frechetMean(mfd, Z)

p = 2
n <- dim(Z)[2]

labels = c("08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18") 
plot(x=NULL,y=NULL,xlim=c(32, 53),ylim=c(32, 38.5), xlab = "longitude", ylab = "latitude", cex.lab=1.4)
obs <- car2sph(x=Z[1,(p+2):n], y=Z[2,(p+2):n], z=Z[3,(p+2):n])
text(obs[,1], obs[,2], labels=labels)
points(obs[,1], obs[,2], cex=3)
lines(obs[,1], obs[,2])


para2 <- getAlpha(Z, mu, model=2, p=2)
n <-dim(para2$U)[2]
pre2 <- predict.hs(Z, para2$U, mu, para2$alpha, model=2, projection=TRUE, PredictorIndex = c((n-p+1):n))

fit2 <- matrix(rep(0, 3*(n-p)), n-p, 3)
for(i in c(1:(n-p))){
  fit2[i,] <- as.vector(predict.hs(Z, para2$U, mu, para2$alpha, model=2, projection=FALSE, PredictorIndex = c(i:(p+i-1)))$predicted)
}

pre = car2sph(x=fit2[,1], y=fit2[,2], z=fit2[,3])
lines(pre[,1], pre[,2], col="blue")
points(pre[,1], pre[,2], col="blue", pch=2, cex = 3)
text(pre[,1], pre[,2], labels=labels, col="blue")


tru <- car2sph(x=truth[1], y=truth[2], z=truth[3])
points(tru[1], tru[2], pch=16, col="black", cex=3)
text(tru[1], tru[2], labels="19", cex=0.9, font=2, col="white")
lines(c(obs[dim(obs)[1],1], tru[1]), c(obs[dim(obs)[1],2], tru[2]))

pred <- as.vector(pre2$predicted)
pred <- car2sph(x=pred[1], y=pred[2], z=pred[3])
points(pred[1], pred[2], pch=17, col="blue", cex = 3)
text(pred[1], pred[2], labels="19", cex=0.9, font=2, col = "white")
lines(c(pre[dim(pre)[1],1], pred[1]), c(pre[dim(pre)[1],2], pred[2]), col="blue")
legend(x = "topleft",          # Position
       legend = c("observed", "fitted", "observed target", "predicted"),  # Legend texts
       pch = c(1, 2, 16, 17),           # Line types
       col = c("black", "blue", "black", "blue")) 






