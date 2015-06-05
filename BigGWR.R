#
# Big GW experiments - in the sandbox
#
# Version 1.0 - Martin Charlton, Maynooth University
# Licensed under GPL V3: https://gnu.org/licenses/gpl.html
#
# Use this at your own risk - this is totally experimental to see
# whether it's possible to run GWR on reasonably large problems
#
# Test system is a Dell Precision T7400 with a 3.16GHz
#   Intel Xeon X5460 CPU running Windows XP SP3. 
#   It has 3.25Gb of RAM
#
# OAC example takes ~400 seconds per iteration for 11 iterations
#   to fit the test model with 3 independent variables and
#   232296 observations
#
# The property price data has 77299 observations and there
#    are 4 independent variabes.  It takes ~140 seconds per
#    iteration for 9 iterations for fit the model
# 
# Issues
# ======
# Kernel: we're restricted to a bisquare kernel at the moment, because
#         I'm using RANN to access a 2-dtree for the nearest neighbour
#         searches
#
# Memory: Because of this we need a matrix of size N x n where n is the
#         upper limit on the number of nearest neighbours in the
#         calibration procedure.  On my XP system this puts a
#         limit of about 250. That's about 58,000,000 cells in the 
#         array.   IF a larger number is required, we'll have to
#         start splitting the problem.  Split into two for 500 nns,
#         four for 1000 nns and so on. 
#
#
library(GWmodel)
library(RANN)
library(rgdal)
library(RColorBrewer) # creates nice color schemes
library(classInt)     # finds class intervals for continuous variables

################################################################################
################################################################################
###
### BIG GWR Functions
###
################################################################################
################################################################################

###
### Modified version of nn2 from RANN
###    - remove distances from data returned
###
### Input:
###   data: coords(spdf) of length N
###   query: coords(spdf)
###   k:     desired number of near neighbours
###
### Output:
###   nn.indexes: matrix of N x k near neighbours - sotred by distance in each row
###
nn3 <- function (data, query = data, k = min(10, nrow(data)), treetype = c("kd",
                        "bd"), searchtype = c("standard", "priority", "radius"),
                        radius = 0, eps = 0)
{
    dimension <- ncol(data)
    ND <- nrow(data)
    NQ <- nrow(query)
    if (ncol(data) != ncol(query))
        stop("Query and data tables must have same dimensions")
    if (k > ND)
        stop("Cannot find more nearest neighbours than there are points")
    searchtypeInt = pmatch(searchtype[1], c("standard", "priority",
        "radius"))
    if (is.na(searchtypeInt))
        stop(paste("Unknown search type", searchtype))
    treetype = match.arg(treetype, c("kd", "bd"))
    if (is.data.frame(data))
        data <- unlist(data, use.names = FALSE)
    if (!is.matrix(query))
        query <- unlist(query, use.names = FALSE)
    results <- .C("get_NN_2Set", as.double(data), as.double(query),
        as.integer(dimension), as.integer(ND), as.integer(NQ),
        as.integer(k), as.double(eps), as.integer(searchtypeInt),
        as.integer(treetype == "bd"), as.double(radius * radius),
        nn.idx = integer(k * NQ), nn = double(k * NQ), PACKAGE = "RANN")
    nn.indexes = matrix(results$nn.idx, ncol = k, byrow = TRUE)
    return(nn.indexes)
}

###
### CV score calibration function
###
### Input:
###   bw:       bandwidth
###   spdf:     spatial data frame containing variables
###   yvar:     index of y variable in spdf@data
###   xvars:    vector of indices of x variables
###
### Output:
###   [cvs]     LOO crossvalidation score
###

cv.test <- function(bw,spdf,yvar,xvars) {
   N <- nrow(spdf)
   coords <- coordinates(spdf)

   if(exists("kd.search")) rm("kd.search")

   kd.search <- nn3(coordinates(spdf),coordinates(spdf),bw)

   betaHat <- matrix(NA,N,length(xvars)+1)

   cv.score <- rep(NA,N)
   for (iter in 1:N) {

      rows.needed <- kd.search[iter,]
      X <- as.matrix(cbind(1,spdf@data[rows.needed,xvars]))
      Y <- spdf@data[rows.needed,yvar]
      #W <- sqrt(rep(1,bw))             # boxcar kernel
      xf <- coords[iter,1]             # bisquare kernel below
      yf <- coords[iter,2]
      xd <- coords[rows.needed,1]
      yd <- coords[rows.needed,2]
      dVec <- ((xf - xd) * (xf - xd)) + ((yf - yd) * (yf - yd))
      W <- (1 - dVec / dVec[bw])       # sqrt(  (1-x^2)^2) is 1-x^2
      W[1] <- 0

      Xw <- sweep(X, 1, W, "*")
      Yw <- matrix(Y * W,bw,1)

      betaHat[iter,] <- solve(t(Xw) %*% Xw, t(Xw) %*% Yw)

      Yhat <-  X[1,] %*% betaHat[iter,]
      cv.score[iter] <- (Y[1] - Yhat)*(Y[1] - Yhat)


      #if(iter %% 1000 == 0) print(iter)
    }
    cvs <- sqrt(sum(cv.score)/N)
    cat(bw, cvs, "\n")
    cvs
}

###
### LOO Calibration by golden section search
###
### Input:
###   spdf:     spatial data frame containing variables
###   yvar:     index of y variable in spdf@data
###   xvars:    vector of indices of x variables
###   xU:       upper limit on bandwdith estimate
###   xL:       lower limit on bandwidth estimate
###
### Output:
###   [xopt]    bandwidth which minimises LOO CV score

calib.big <- function(spdf,yvar,xvars,xU,xL) {
    ptm <- proc.time()[3]
    adapt.bw <- TRUE
    eps = 1e-04
    R <- (sqrt(5) - 1)/2
    iter <- 1
    d <- R * (xU - xL)
    if (adapt.bw) {
        x1 <- floor(xL + d)
        x2 <- round(xU - d)
    } else {
        x1 <- xL + d
        x2 <- xU - d
    }
    f1 <- cv.test(x1,spdf,yvar,xvars)
    f2 <- cv.test(x2,spdf,yvar,xvars)
    d1 <- f2 - f1
    if (f1 < f2) {
        xopt <- x1
    } else {
        xopt <- x2
    }
    ea <- 100
    while ((abs(d) > eps) && (abs(d1) > eps)) {
        d <- R * d
        if (f1 < f2) {
            xL <- x2
            x2 <- x1
            if (adapt.bw) {
                x1 <- round(xL + d)
            } else {
            x1 <- xL + d
            }
            f2 <- f1
            f1 <- cv.test(x1,spdf,yvar,xvars)
        } else {
            xU <- x1
            x1 <- x2
            if (adapt.bw) {
                x2 <- floor(xU - d)
            } else {
             x2 <- xU - d
            }
            f1 <- f2
            f2 <- cv.test(x2,spdf,yvar,xvars)
        }
        iter <- iter + 1
        if (f1 < f2) {
            xopt <- x1
        } else {
            xopt <- x2
        }
        d1 <- f2 - f1
    }
    ptm <- proc.time()[3] - ptm
    cat("Elapsed time",ptm,"seconds\n")
    xopt
}

###
### GWR fitting function
###
### Input:
###   spdf:     spatial data frame containing variables
###   yvar:     index of y variable in spdf@data
###   xvars:    vector of indices of x variables
###   bw:       bandwidth
###
### Output:
###   [betahat] N x p matrix of parameter estimates
###

gwr.big <- function(spdf,yvar,xvars,bw) {

   N <- nrow(spdf)
   coords <- coordinates(spdf)

   if(exists("kd.search")) rm("kd.search")

   kd.search <- nn3(coordinates(spdf),coordinates(spdf),bw)

   betaHat <- matrix(NA,N,length(xvars)+1)

   for (iter in 1:N) {

      rows.needed <- kd.search[iter,]
      X <- as.matrix(cbind(1,spdf@data[rows.needed,xvars]))
      Y <- spdf@data[rows.needed,yvar]
      #W <- sqrt(rep(1,bw))             # boxcar kernel
      xf <- coords[iter,1]             # bisquare kernel below
      yf <- coords[iter,2]
      xd <- coords[rows.needed,1]
      yd <- coords[rows.needed,2]
      dVec <- ((xf - xd) * (xf - xd)) + ((yf - yd) * (yf - yd))
      W <- (1 - dVec / dVec[bw])       # sqrt(  (1-x^2)^2) is 1-x^2

      Xw <- sweep(X, 1, W, "*")
      Yw <- matrix(Y * W,bw,1)

      betaHat[iter,] <- solve(t(Xw) %*% Xw, t(Xw) %*% Yw)
   }

      betaHat
}




########################################################################
########################################################################
###
### First experiment - UK census data
###
########################################################################
########################################################################

setwd("H:/OAC_2011")    # office
setwd("E:/OAC_2011")    # laptop


PopLkp <- read.csv("2011 OA Population and Lookup.csv",stringsAsFactors=F)

#> dim(PopLkp)
#[1] 232296      8
#> head(PopLkp)
#         OA LATITUDE  LONGITUDE LOCAL_AUTHORITY_CODE LOCAL_AUTHORITY_NAME
#1 E00000001 51.52027 -0.0950106            E09000001       City of London
#2 E00000003 51.51989 -0.0968518            E09000001       City of London
#3 E00000005 51.51903 -0.0962386            E09000001       City of London
#4 E00000007 51.51686 -0.0982424            E09000001       City of London
#5 E00000010 51.52240 -0.0973724            E09000001       City of London
#6 E00000012 51.52209 -0.0960576            E09000001       City of London
#
#  REGION_OR_COUNTRY_CODE REGION_OR_COUNTRY_NAME POPULATION
#1              E12000007                 London        194
#2              E12000007                 London        250
#3              E12000007                 London        367
#4              E12000007                 London        123
#5              E12000007                 London        102
#6

OACv60 <- read.csv("2011 OAC 60 variables/2011_OAC_Raw_kvariables.csv",stringsAsFactors=F)
dim(OACv60)
#OACvars60lkp <- read.csv("2011 OAC 60 variables/2011_OAC_Raw_kvariables_Lookup.csv",stringsAsFactors=F)
#dim(OACvars60lkp)

Educ   <- OACv60$k039
Unemp  <- OACv60$k045
SocHou <- OACv60$k032
PubTra <- OACv60$k042

rm(OACv60)


OACspdfGeo <- SpatialPointsDataFrame(data.frame(PopLkp$LONGITUDE,PopLkp$LATITUDE),data.frame(PopLkp$OA,Educ,Unemp,SocHou,PubTra))
proj4string(OACspdfGeo) <- CRS("+init=epsg:4326")
OACspdf <- spTransform(OACspdfGeo,CRS("+init=epsg:27700"))
plot(OACspdf,pch=16,cex=0.25,col="red")
rm(OACspdfGeo)

cor(OACspdf@data[,2:5])

plotvar <- OACspdf$Educ
nclr <- 8
plotclr <- brewer.pal(nclr,"BuPu")
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
plot(OACspdf,pch=16,cex=0.25,col=colcode)
title("Level 4 Qual 8-Quantiles")






###
### Calibration results
###

    yvar <- 2
    xvars <- 3:5
    xU <- 240
    xL <- 30
    bw.opt <- calib.big(OACspdf,yvar,xvars,xL,xU)
    print(bw.opt)



bw.vec <- c(     159,      110,       78,       60,       47,       66,       72,       64,       68,
                  65,       67,       65,       66)
cv.vec <- c(21.83061, 21.53984, 21.39067, 21.38754, 21.53422, 21.36987, 21.37529, 21.37240, 21.37030,
            21.37037, 21.37131, 21.37037, 21.36987)

N.n    <- 232296
bw.opt <- 66

betaHat <- gwr.big(OACspdf,yvar,xvars,bw.opt)

ord.vec <- order(bw.vec)
cvc     <- splinefun(bw.vec[ord.vec], cv.vec[ord.vec])
x.cvc   <- seq(min(bw.vec),max(bw.vec),l=100)
y.cvc   <- cvc(x.cvc)
plot(x.cvc,y.cvc,type="l",xlab="Bandwidth",ylab="Crossvalidation Score")
title("Calibration for 'Big Data' GWR: N=232296")
points(bw.vec[ord.vec], cv.vec[ord.vec],pch=16,col="red",cex=0.5)




##
## Fit 'optimal' model
##

plotvar <- betaHat[,4]
nclr <- 8
plotclr <- (brewer.pal(nclr,"BuPu"))   # rev() not for 4 or 1
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
plot(OACspdf,pch=16,cex=0.25,col=colcode)
title("Level 4 Qual Public Transport Coefficient")
legend(-30000, 1250000, legend=names(attr(colcode, "table")),
    fill=attr(colcode, "palette"), cex=0.6, bty="n")





##############################################################################
##############################################################################
###
### Experiment two - housing data
###
##############################################################################
##############################################################################

setwd("H:/OAC_2011")    # office
Data90 <- read.table("GWRmerge90.dat",header=F,sep=" ",stringsAsFactors=F)
colnames(Data90) <- c("Easting","Northing","Purprice","BldIntWr","BldPostW","Bld60s",
 "Bld70s","Bld80s","TypDetch","TypSemiD","TypFlat","GarSingl","GarDoubl","Tenfree",
 "CenHeat","BathTwo","BedTwo","BedThree","BedFour","BedFive","NewPropD","FlorArea",
 "NoCarHh","CarspP","ProfPct","UnskPct","RetiPct","Saleunem","Unemploy","PopnDnsy")

D90spdf <- SpatialPointsDataFrame(data.frame(Data90$Easting, Data90$Northing),
                                  data.frame(Data90$Easting, Data90$Northing,Data90$Purprice,
                                             Data90$FlorArea,Data90$TypDetch,Data90$TypSemiD,Data90$TypFlat))

yvar  <- 3
xvars <- 4:7
xL    <- 250
xU    <- 500

bw.opt  <- calib.big(D90spdf,yvar,xvars,xL,xU)

betaHat <- gwr.big(D90spdf,yvar,xvars,bw.opt)

###
### Plot some parameter estimates
###
data(EWOutline)
plotvar <- betaHat[,2]
nclr    <- 8
plotclr <- rev(brewer.pal(nclr,"Spectral"))
class   <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
plot(D90spdf,pch=16,cex=0.25,col=colcode)
plot(ewoutline,add=T)
title("NW Housing Data: Floor Area Parameter Variation")

plot(Data90$FlorArea,Data90$Purprice,pch=16,cex=0.25,col=colcode)

