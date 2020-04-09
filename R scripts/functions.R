
# get the distance between two points on the globe.
#
# args:
# lat1 - latitude of first point.
# long1 - longitude of first point.
# lat2 - latitude of first point.
# long2 - longitude of first point.
# radius - average radius of the earth in km
#
# see: http://en.wikipedia.org/wiki/Great_circle_distance

greatCircleDistance <- function(lat1, long1, lat2, long2, radius=6372.795){
  sf <- pi/180
  lat1 <- lat1*sf
  lat2 <- lat2*sf
  long1 <- long1*sf
  long2 <- long2*sf
  lod <- abs(long1-long2)
  radius * atan2(
    sqrt((cos(lat1)*sin(lod))**2 +
           (cos(lat2)*sin(lat1)-sin(lat2)*cos(lat1)*cos(lod))**2),
    sin(lat2)*sin(lat1)+cos(lat2)*cos(lat1)*cos(lod)
  )
}

# Calculate the nearest point using latitude and longitude.
# and attach the other args and nearest distance from the
# other data.frame. 

dist.merge <- function(x, y, xlongnme, xlatnme, ylongnme, ylatnme){
  tmp <- t(apply(x[,c(xlongnme, xlatnme)], 1, function(x, y){
    dists <- apply(y, 1, function(x, y) greatCircleDistance(x[2],
                                                            x[1], y[2], y[1]), x)
    cbind(1:nrow(y), dists)[dists == min(dists),,drop=F][1,]
  }
  , y[,c(ylongnme, ylatnme)]))
  tmp <- cbind(x, min.dist=tmp[,2], y[tmp[,1],-match(c(ylongnme,
                                                       ylatnme), names(y))])
  row.names(tmp) <- NULL
  tmp
}


#panel correlation plots to assess data distribution
panel.hist <- function(x, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(usr[1:2], 0, 1.5) )     
  h <- hist(x, plot = FALSE)     
  breaks <- h$breaks; nB <- length(breaks)     
  y <- h$counts; y <- y/max(y)     
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...) 
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(0, 1, 0, 1))     
  r <- abs(cor(x, y, use = "complete"))   
  txt <- format(c(r, 0.123456789), digits = digits)[1]     
  txt <- paste0(prefix, txt)     
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)     
  text(0.5, 0.5, txt, cex = cex.cor * r) }

#This code is modified for D.L. Miller's dsm package for distance sampling, from
#the rqgam.check function. The code is designed to extract randomized quantile 
#residuals from GAMs, using the family definitions in mgcv. Note statmod only
#supports RQ residuals for the following families: Tweedie, Poisson, Gaussian,  Any errors are due to Eric Pedersen
library(statmod) #This has functions for randomized quantile residuals
rqresiduals = function (gam.obj) {
  if(!"gam" %in% attr(gam.obj,"class")){
    stop('"gam.obj has to be of class "gam"')
  }
  if (!grepl("^Tweedie|^Negative Binomial|^poisson|^binomial|^gaussian|^Gamma|^inverse.gaussian",
             gam.obj$family$family)){
    stop(paste("family " , gam.obj$family$family, 
               " is not currently supported by the statmod library, 
                 and any randomized quantile residuals would be inaccurate.",
               sep=""))
  }
  if (grepl("^Tweedie", gam.obj$family$family)) {
    if (is.null(environment(gam.obj$family$variance)$p)) {
      p.val <- gam.obj$family$getTheta(TRUE)
      environment(gam.obj$family$variance)$p <- p.val
    }
    qres <- qres.tweedie(gam.obj)
  }
  else if (grepl("^Negative Binomial", gam.obj$family$family)) {
    if ("extended.family" %in% class(gam.obj$family)) {
      gam.obj$theta <- gam.obj$family$getTheta(TRUE)
    }
    else {
      gam.obj$theta <- gam.obj$family$getTheta()
    }
    qres <- qres.nbinom(gam.obj)
  }
  else {
    qres <- qresid(gam.obj)
  }
  return(qres)
}
