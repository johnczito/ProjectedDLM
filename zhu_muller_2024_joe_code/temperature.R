library(MASS)
library(readxl)
library(pracma)
library(fdapace)
library(manifold)

### JFK 
indexes <- read.csv("jfk.csv")
indexes <- as.data.frame(indexes)

date <- indexes$DATE
DateArray <- sapply(date, function(s){
  s<-as.character(s)
  return(strtoi(unlist(strsplit(s, "-")), 10L))  
})

year <- DateArray[1,]
month <- DateArray[2,]
index <- which(month %in% c( 6, 7, 8, 9))
indexes <- indexes[index, ]
year <- year[index]
yr <- seq(1990, 2019, by=1)

dens <- vector(mode = "list", length = length(yr))
l1 <- floor(quantile(indexes$TMIN, 0.01, na.rm = TRUE))
l2 <- floor(quantile(indexes$TMIN, 0.99, na.rm = TRUE))
u1 <- floor(quantile(indexes$TMAX, 0.01, na.rm = TRUE))
u2 <- floor(quantile(indexes$TMAX, 0.99, na.rm = TRUE))

xgrid <- c(l1:l2)
ygrid <- c(u1:u2)
xtemp <- rep(xgrid, length(ygrid))
ytemp <- as.vector(sapply(ygrid, function(z){
  return(rep(z, length(xgrid)))
}))
xin <- cbind(xtemp, ytemp)
N <- length(xtemp)
xout1 = seq(l1, l2, length.out=51)
xout2 = seq(u1, u2, length.out=51)
h1 <- diff(range(xgrid))/5
h2 <- diff(range(ygrid))/5
a <- diff(range(xout1))/50*diff(range(xout2))/50

for(t in c(1:length(yr))){
  #index <- which(year %in% c(yr[t]:(yr[t]+2)))
  index <- which(year %in% c(yr[t]))
  temp <- indexes[index, ]
  
  x <- temp$TMIN
  y <- temp$TMAX
  
  index <- which(!is.na(x))
  x <- x[index]
  y <- y[index]
  
  index <- which(!is.na(y))
  x <- x[index]
  y <- y[index]
  
  yin <- rep(0, N)
  for (i in c(1:N)) {
    index1 <- which(x==xin[i,1])
    index2 <- which(y==xin[i,2])
    yin[i] <- length(intersect(index1, index2))
  }
  #dens[[t]] <- kde2d(x, y, h =c(h1, h2), n=51, lims = c(l1, l2, u1, u2))

  z <- Lwls2D(bw=c(h1, h2),
    kern = "epan",
    xin=xin,
    yin=yin,
    win = NULL,
    xout1 = xout1,
    xout2 = xout2,
    xout = NULL,
    subset = NULL,
    crosscov = TRUE
  )
  index <- which(z<0)
  z[index] <-0
  z <- z/sum(z)
  z <- z/a
  dens[[t]] <- list(x=xout1, y=xout2, z=z)
}

names = c(1:30)-1
for(i in 1:30){
  mypath <- file.path("figure2d/",paste("2d-", names[i], ".png", sep = ""))
  png(file=mypath)
  mytitle = paste("Year", yr[i])
  dat <- list("x" = dens[[i]]$x, "y"=dens[[i]]$y, "z"=dens[[i]]$z)
  filled.contour3(dat, xlim = c(l1, l2), ylim = c(u1, u2), zlim = c(0, 0.007), cex.axis=1.5, cex.lab=1.5, plot.title =title(main=mytitle))
  dev.off()
}


L <- length(yr)
Z <- sapply(c(1:(L-1)), function(i){
  temp <- dens[[i]]
  temp$z <- temp$z*a
  #temp$z <- temp$z/sum(temp$z)
  z <- sqrt(as.vector(temp$z))
  return(z)
})

zmax <- quantile(sapply(c(1:length(dens)), function(i){
  temp <- dens[[i]]
  return(max(temp$z))
}), 0.99)

mfd <- structure(1, class='Sphere')
mu <- frechetMean(mfd, Z)

# model 1
p=5
para1 <- getAlpha(Z, mu, model=1, p=p)
n <-dim(para1$U)[2]
pre1 <- predict.hs(Z, para1$U, mu, para1$alpha, model=1, projection=TRUE,  PredictorIndex = c((n-p+1):n))
dat <- list("x" = dens[[1]]$x, "y"=dens[[1]]$y, "z"=matrix(pre1$predicted^2/a, 51, 51))
filled.contour3(dat, xlim = c(l1, l2), ylim = c(u1, u2), zlim = c(0, zmax), cex.axis=1.5, cex.lab=1.5)

# model 2
p=5
para2 <- getAlpha(Z, mu, model=2, p=p)
n <-dim(para2$U)[2]
pre2 <- predict.hs(Z, para2$U, mu, para2$alpha, model=2, projection=TRUE, PredictorIndex = c((n-p+1):n))
dat <- list("x" = dens[[1]]$x, "y"=dens[[1]]$y, "z"=matrix(pre2$predicted^2/a, 51, 51))
filled.contour3(dat, xlim = c(l1, l2), ylim = c(u1, u2), zlim = c(0, zmax), cex.axis=1.5, cex.lab=1.5)
filled.contour3(dens[[length(dens)]], xlim = c(l1, l2), ylim = c(u1, u2), zlim = c(0, zmax), cex.axis=1.5,cex.lab=1.5)

