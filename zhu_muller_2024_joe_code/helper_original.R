
cv.hs <- function(Z, model){
  d <- dim(Z)[1]
  n <- dim(Z)[2]
  
  n.test <- floor(n/10)
  n.train <- n - n.test
  
  ERR <- sapply(c(1:n.test), function(t){
    truth <- Z[,c(n-t+1)]
    Z.train <- Z[,c((n-t-n.train+1):(n-t))]
    
    
      
    mfd <- structure(1, class='Sphere')
    mu <- tryCatch(RFPCA:::frechetMean(mfd, Z), error=function(e) return(rowMeans(Z.train)/sqrt(sum((rowMeans(Z.train))^2))))
    
    err <- sapply(c(1:5), function(p){
      para <- getAlpha(Z.train, mu, model=model, p=p)
      n.temp <-dim(para$U)[2]
      pre <- predict.hs(Z.train, para$U, mu, para$alpha, model=model, projection=FALSE,  PredictorIndex = c((n.temp-p+1):n.temp))
      return(acos(sum(pre$predicted*truth)))
    })
    
    return(err)
  })
  
  return(which.min(rowMeans(ERR)))
} 

getAlpha <- function(Z, mu, model, p){
  
  U <- getU(Z, mu, model)
  
  d <- dim(Z)[1]
  n <- dim(U)[2]
  
  V <- matrix(rep(0, d*n), d, n)
  if(model==1){
    V <- matrix(rep(mu, n), ncol = n)
  }else{
    V<- Z[, c(1:n)]
  }
  
  D <- matrix(rep(0, n*n), n, n)
  for (i in c(1:n)) {
    for (j in c(1:n)) {
      D[i,j] <- 2*sum(V[,i]*V[,j])*sum(U[,i]*U[,j])-2*sum(V[,i]*U[,j])*sum(V[,j]*U[,i])
    }
  }
  M <- matrix(rep(rowMeans(D), n), ncol = n)
  D <- D - M - t(M) + mean(D)
  
  sig <- rep(0, p)
  for (i in c(1:p)) {
    sig[i]<- mean(diag(D[c((i+1):n), c(1:(n-i))]))
  }
  sig0 <- mean(diag(D))
  
  A <- matrix(rep(0, p*p), p, p)
  for (i in c(1:p)) {
    for (j in c(1:p)) {
      if(i==j){
        A[i,j] <- sig0
      }else{
        A[i,j] <- sig[abs(i-j)]
      }
    }
  }
  
  alpha <- solve(A)%*%sig
  
  par <- list("Z" = Z, "U"=U, "alpha"=alpha, "model"=model)
  return(par)
}

getU <- function(Z, mu=NULL, model){
  n <- dim(Z)[2]
  
  if(model==1){
    U <- sapply(c(1:n), function(i){
      u1 <- mu
      u2 <- Z[,i] -  sum(Z[,i]*u1)*u1
      u2 <- u2/sqrt(sum(u2^2))
      theta <- acos(sum(mu*Z[,i]))
      return(theta*u2)
    })
  }else{
    U <- sapply(c(1:(n-1)), function(i){
      u1 <- Z[,i]
      u2 <- Z[,i+1] -  sum(Z[,i+1]*u1)*u1
      u2 <- u2/sqrt(sum(u2^2))
      theta <- acos(sum(Z[,i]*Z[,i+1]))
      return(theta*u2)
    })
  }
}

applyS <- function(v, Z, U, mu, alpha, model, PredictorIndex){
  p <- length(alpha)
  if(model==1){
    R <- sapply(c(1:dim(U)[2]), function(i){
      return(sum(mu*v)*U[,i] - sum(U[,i]*v)*mu)
    }) 
  }else{
    R <- sapply(c(1:dim(U)[2]), function(i){
      return(sum(Z[,i]*v)*U[,i] - sum(U[,i]*v)*Z[, i])
    })
  }
  
  B0 <- rowMeans(R)
  b <- B0*(1-sum(alpha))
  if(p>1){
    a <- R[, PredictorIndex]%*%alpha[c(p:1)]
  }else{
    a <- R[, PredictorIndex]*alpha[c(p:1)]
  }
  
  return(b+a)
}

predict.hs <- function(Z, U, mu, alpha, model, projection, PredictorIndex){
  p <- length(alpha)
  
  if(model==1){
    v=mu
  }else{
    v <- Z[, (PredictorIndex[p]+1)]
  }
  
  d <- length(v)
  N=10
  V <- matrix(rep(0, d*N), ncol = N)
  V[,1] <- v
  for(i in c(2:N)){
    V[,i] <- applyS(V[,(i-1)], Z, U, mu, alpha, model, PredictorIndex)
  }
  
  coef <- sapply(c(1:N), function(i){
    return(1/factorial(i-1))
  })
  
  re <- V%*%coef
  
  ##########
  # diag(3)+S + S%*%S/factorial(2) + S%*%S%*%S/factorial(3) +  S%*%S%*%S%*%S/factorial(4) +
  # S%*%S%*%S%*%S%*%S/factorial(5) + S%*%S%*%S%*%S%*%S%*%S/factorial(6) + S%*%S%*%S%*%S%*%S%*%S%*%S/factorial(7)
  
  #### project to the positive quatant
  # if(sum(re<0)>0){
  #   index <- which(re<0)
  #   re[index] = 0
  #   re <- re/sqrt(sum(re^2))
  # }
  
  
  re <- re/sqrt(sum(re^2))
  ra <- 1
  # if(projection==TRUE){
  #   if(sum(re<0)>0){
  #     u1 <- v
  #     u2 <- re -  sum(re*u1)*u1
  #     u2 <- u2/sqrt(sum(u2^2))
  #     theta <- acos(sum(v*re))
  #     u2 <- u2*theta
  #     right <- ra
  #     left=0
  #     
  #     V[,1] <- v
  #     for(i in c(2:N)){
  #       V[,i] <- sum(u1*V[,i-1])*u2 - sum(u2*V[,i-1])*u1
  #     }
  #     
  #     K= 10
  #     while (K>1) {
  #       
  #       ra <- (left+right)/2
  #       coef <- sapply(c(1:10), function(i){
  #         return((ra)^(i-1)/factorial(i-1))
  #       })
  #       re <- V%*%coef
  #       
  #       if(sum(re<0)>0){
  #         right=ra
  #       }else{
  #         left=ra
  #       }
  #       K=K-1
  #     }
  #     ra <- left
  #     coef <- sapply(c(1:N), function(i){
  #       return((ra)^(i-1)/factorial(i-1))
  #     })
  #     re <- V%*%coef
  #   }
  # }
  return(list("predicted"=re, "ratio"=ra))
}

proj<- function(re, v){
  ra=1
  if(sum(re<0)>0){
    u1 <- v
    u2 <- re -  sum(re*u1)*u1
    u2 <- u2/sqrt(sum(u2^2))
    theta <- acos(sum(v*re))
    u2 <- u2*theta
    right <- ra
    left=0
    
    N=10
    d <- length(v)
    V <- matrix(rep(0, d*N), ncol = N)
    V[,1] <- v
    for(i in c(2:N)){
      V[,i] <- sum(u1*V[,i-1])*u2 - sum(u2*V[,i-1])*u1
    }

    K= 10
    while (K>1) {
      ra <- (left+right)/2
      coef <- sapply(c(1:10), function(i){
        return((ra)^(i-1)/factorial(i-1))
      })
      re <- V%*%coef

      if(sum(re<0)>0){
        right=ra
      }else{
        left=ra
      }
      K=K-1
    }
    ra <- left
    coef <- sapply(c(1:N), function(i){
      return((ra)^(i-1)/factorial(i-1))
    })
    re <- V%*%coef
  }
  return(re)
}

density2d <- function(x1, x2, xout1=NULL, xout2=NULL){
  ## xout1 and xout2 are assumed to be equi-distanced
  if(is.null(xout1) || is.null(xout2)){
    l1 <- min(x1)
    l2 <- max(x1)
    u1 <- min(x2)
    u2 <- max(x2)
    
    xout1 = seq(l1, l2, length.out=51)
    xout2 = seq(u1, u2, length.out=51)
  }
  h1 <- diff(xout1)[1]
  h2 <- diff(xout2)[1]
  a <- h1*h2
  
  xtemp <- rep(xout1, length(xout2))
  ytemp <- as.vector(sapply(xout2, function(z){
    return(rep(z, length(xout1)))
  }))
  xin <- cbind(xtemp, ytemp)
  N <- dim(xin)[1]
  
  index <- which(!is.na(x1))
  x1 <- x1[index]
  x2 <- x2[index]
  
  index <- which(!is.na(x2))
  x1 <- x1[index]
  x2 <- x2[index]
  
  yin <- rep(0, N)
  for (i in c(1:N)) {
    index1 <- which((x1 >= (xin[i,1] - h1/2)) & (x1 < (xin[i,1] + h1/2)))
    index2 <- which((x2 >= (xin[i,2] - h2/2)) & (x2< (xin[i,2] + h2/2)))
    yin[i] <- length(intersect(index1, index2))
  }
  
  k1 <- diff(range(xout1))/5
  k2 <- diff(range(xout2))/5
  z <- Lwls2D(bw=c(k1, k2),
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
  
  return(list(x=xout1, y=xout2, z=z))
}

filled.contour3 <-function (x = seq(0, 1, length.out = nrow(z)),
                            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page
  
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  # on.exit(par(par.orig))
  # w <- (3 + mar.orig[2]) * par("csi") * 2.54
  # par(las = las)
  # mar <- mar.orig
  plot.new()
  # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                  col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

filled.legend <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                                        length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
                           ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                           levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                           col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                           key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                           axes = TRUE, frame.plot = axes, ...) 
{
  # modification of filled.contour by Carey McGilliard and Bridget Ferris
  # designed to just plot the legend
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #  on.exit(par(par.orig))
  #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  #  par(las = las)
  #  mar <- mar.orig
  #  mar[4L] <- mar[2L]
  #  mar[2L] <- 1
  #  par(mar = mar)
  # plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  
  if (frame.plot)
    box()
  if (missing(plot.title))
    title(...)
  else plot.title
  invisible()
  box()
}
