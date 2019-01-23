hummingage <- function(path,infile,segments,DiscardOutliers,EstimateNullRate,token,plot_format) {
    ##
    ## usage:
    ## hummingage("/home/user/data/","myfile.csv",c(NaN),TRUE,"estimate","KJmpmGgQ","png,pdf")
    ## hummingage("/home/user/data/","myfile.csv",c(124),FALSE,"estimate","KJmpmGgQ","png")
    ## hummingage("/home/user/data/","myfile.csv",c(334,856,1223),TRUE,"estimate","KJmpmGgQ","pdf")
    ## segments in cm

    ## path: (string): is the path to the infile
    ## infile: (string): is the data file according to a certain format
    ## segments: (integer array): is an array indicating roughly the positions of the segments, use c(NaN) for no segments
    ## DiscardOutliers: (boolean): can be TRUE or FALSE. If true, a simple algorithm is applied to identify outliers
    ## EstimateNullRate: (string or integer): possible values are
    ##                   "estimate": the x-axis offset is estimated from the data
    ##                   "null": the x-axis offset is set to 0
    ##                   negative or positive integer: x-axis offset is set to the inserted value, e.g. -42, 14, ... (in cm)
    ## token: (string): this is a random string, which can be choosen arbitrarily
    ## plot_format (string): this can be "png,pdf", "png" or "pdf" and defines what kinds of images are written to disk
    ##                       both, png and pdf or only png or pdf
    ##


    ## full file-path
    xfile <- paste(path,infile,sep="");

    ## transform plot_format to array
    plot_format_tmp <- strsplit(plot_format,",")
    plot_format <- plot_format_tmp[[1]]

    
    ## load the infile and save as ageDataOrig
    ageDataOrig <- read.table(xfile, header = TRUE,sep=",")

    ## if no segment will be used set it to empty
    if (is.nan(segments)){
        segments = c()
    }

    ## calculate length and estimate rates
    lengthData <- length(ageDataOrig$depth)
    diffDepths <- diff(ageDataOrig$depth)
    accRatesOrig <- diff(ageDataOrig$age)/diffDepths

    ##--------------------------------------------------------------------------------------------##
    ## this block is only executed if DiscardOutliers == TRUE
    ## identify negative rates, i.e. outliers
    ## the difficulty is to identify the outlier
    ## i.e. which age is the outlier, because the rate is calculated from two ages
    ## thus we take the ages left and right to the two ages which make the negative rate
    ## and then we calculate which age is further away from the mean of the two left
    ## and rights
    if (DiscardOutliers==TRUE){
        ## find negative rates
        negRatesIndex <- which(accRatesOrig < 0)
    } else {
        ## set this to empty
        negRatesIndex <- integer(0)
    }
    ## make a copy of the indices of the negative rates
    negRatesIndexCopy <- negRatesIndex

    ## go through the negative rates and decide which age is the outlier, left or right
    m <- 1
    for (i in negRatesIndexCopy) {
        ## skip the first
        if (i==1){
            next
        }
        ## calculate the mean between the left and right age next to the two
        ## ages, which where used to calculate the rate
        TheMean <- (ageDataOrig$age[(i-1)]+ageDataOrig$age[(i+2)])/2
        ## distance of the left age to the mean
        dist1 <- abs(ageDataOrig$age[(i)]-TheMean)
        ## distance of the right age to the mean
        dist2 <- abs(ageDataOrig$age[(i+1)]-TheMean)
        ## if the right is further away then this is the outlier
        ## so increment the index
        ## otherwise it's the left and leave the index
        if (dist2>dist1){
            negRatesIndex[m] <- negRatesIndex[m]+1
        }
        m <- m+1
    }
    ##--------------------------------------------------------------------------------------------##

    ## remove extreme rates and exclude the first
    ## because the first must be estimated later below
    ## remove outliers from original data
    ## or copy data as they are, if there are no outliers
    LnegRatesIndex <- length(negRatesIndex)
    if (LnegRatesIndex>0){
        ## remove data
        ageData <- ageDataOrig[-negRatesIndex,]
    } else {
        ## leave data as they are
        ageData <- ageDataOrig
    }

    ## now, outliers have been removed calculate again length, rates etc.
    lengthData <- length(ageData$depth)
    diffDepths <- diff(ageData$depth)
    accRates <- diff(ageData$age)/diffDepths
    sigmaRates <- sqrt((ageData$error[2:lengthData]^2+ageData$error[1:(lengthData-1)]^2))/diffDepths

    ## convert segments from depth to point number
    ## find the position, i.e. between the two closest points
    Lsegments <- length(segments)
    if (Lsegments == 0) {
        segPos = c()
    } else {
        segPos <- rep(0,Lsegments)
        for (i in seq(1,Lsegments)){
            ## difference between depth variable and a sequence position (in cm)
            Xdiff <- ageData$depth-segments[i]
            ## find the indices smaller than the segment position
            a <- which(Xdiff<0)
            ## set the position at the last negative
            segPos[i] <- max(a)
        } 
    }

    ## copy of segPos
    segPosPlot <- segPos
    ## add left and right boundaries to segPos
    segPos <- c(0,segPos,lengthData)
    segPosLength <- length(segPos)

    ##--------------------------------------------------------------------------------------------##
    ##  estimate the x-axis crossing to estimate the first rate
    ##  if we would take the origin 0,0 than the first rate will be wrong
    ##  because there is an offset
    ##
    ##  so, we fit the first segment and estimate the x-crossing
    ##  cut out first segment
    xcrossing <- 0
    if (EstimateNullRate == "estimate"){
        N <- 3
        Lu <- length(1:N)
        u <- ageData$age[1:N]
        X <- replicate(2,rep(1,Lu))
        X[,2] <- ageData$depth[1:N]
        ##  print(X)
        beta <- solve(t(X) %*% X) %*% t(X) %*% u
        ##  print(beta)
        ##  only if slope is positive
        if (beta[2]>0){
            xcrossing <- beta[1]/beta[2]*(-1)
            ##  print(xcrossing)
        } else {
            xcrossing <- 0
        }
    }
    if (EstimateNullRate == "null"){
        xcrossing <- 0
    }
    if (is.numeric(EstimateNullRate) == TRUE){
        ##  then this is xcrossing
        xcrossing <- EstimateNullRate
    }
    NullDiff <- ageData$depth[1]-xcrossing
    NullRate <- ageData$age[1]/NullDiff
    NullError <- ageData$error[1]/NullDiff
    ##--------------------------------------------------------------------------------------------##

    ## now we have the first rate and add it to the rates
    accRates <- c(NullRate,accRates)
    ## print(accRates)
    ## add also a std for the first
    sigmaRates <- c(NullError,sigmaRates)

    ## initialise arrays for the fits
    fittedRates <- rep(0,lengthData)
    fittedError <- rep(0,lengthData)

    ##--------------------------------------------------------------------------------------------##
    ## fitting the segments
    for (i in seq(2,segPosLength)) {
        ##  least square
        Xsegment <- (segPos[i-1]+1):segPos[i]
        Lsegment <- length(Xsegment)
        ##  create X matrix
        X <- replicate(2,rep(1,Lsegment))
        X[,2] <- ageData$depth[Xsegment]
        ##  create error matrix
        g <- solve(t(X) %*% X)
        ##  print(g)
        beta <- g %*% t(X) %*% accRates[Xsegment]
        ##  print(beta)
        fittedRates[Xsegment] <- X %*% beta
        ##  res
        fittedError[Xsegment] <- sqrt(1/(Lsegment-1) * sum((fittedRates[Xsegment]-accRates[Xsegment])^2))
        ##  print(fittedError[Xsegment])
    }
    ##--------------------------------------------------------------------------------------------##

    
    ##--------------------------------------------------------------------------------------------##
    ##  BIC       
    ##  BIC for Gaussian approx. is BIC=N*ln(epsilon^2)+k*ln(N)
    ##  for fitting a linear regression one need k=q+2 parameters to be estimated
    ##  one for the intercept, q for the slopes and one for the rmse
    ##  thus in the case of fitting the mean I need k=2 (intercept, rmse)
    ##  BIC[segments] <- N * log(rmse^2) + 2*segments * log(N)
    mse <- 1/lengthData * sum((fittedRates-accRates)^2)
    BIC <- lengthData * log(mse) + (2*(segPosLength-1)+1) * log(lengthData)
    ## (segPosLength-1) because above I have added 0,segPos,N to segPos
    ## thus (segPosLength-1) is exactly the number of segments
    ## write out the BIC
    write(BIC,file=paste(path,token,"_BIC.txt",sep=""),ncolumns=1)
    ##--------------------------------------------------------------------------------------------##


    ##--------------------------------------------------------------------------------------------##
    ## apply the Bayesian approach, i.e. for the Gaussian assumption the weighted mean
    BayesFittedRates <- (fittedRates * sigmaRates^2 + accRates * fittedError^2)/(fittedError^2 + sigmaRates^2)
    ## disable Bayes
    ##BayesFittedRates <- fittedRates
    
    ## and the weighted std
    BayesSigmaRates <- sqrt((fittedError^2 * sigmaRates^2)/(fittedError^2 + sigmaRates^2))
    ## disable Bayes
    ##BayesSigmaRates <- fittedError

    
    ## calculate the ages from the rates
    ageFit <- cumsum(BayesFittedRates * c(NullDiff,diffDepths))
    ##  print("sum1")
    ##  print(sum((ageFit-ageData$age)^2))

    ##  and the errors
    ageFitError <- rep(0,lengthData)
    ageFitError[1] <- BayesSigmaRates[1]*NullDiff
    for (i in seq(1,(lengthData-1))) {
        ageFitError[i+1] <- sqrt(ageFitError[i]^2+(BayesSigmaRates[i+1]*diffDepths[i])^2)
    }
    ##--------------------------------------------------------------------------------------------##
    

    ##--------------------------------------------------------------------------------------------##
    ## plotting the intermediate result of fitted segments
    for (i in plot_format){
        if (i == "png"){
            png(filename=paste(path,token,"_intermediate.png",sep=""),width=2000,height=1400,pointsize=10)#,width=11,height=8)
        }
        if (i == "pdf"){
            pdf(file=paste(path,token,"_intermediate.pdf",sep=""),width=40,height=28)#,width=11,height=8)
        }
        
        par(mar=c(14,14,4,2),mgp=c(10,4,0),bg=rgb(248/255,250/255,252/255))

        
        ymin <- min(c(accRates-sigmaRates),na.rm=TRUE)
        ymax <- max(c(accRates+sigmaRates),na.rm=TRUE)
        xmin <- 0 ##min(ageDataOrig$depth)
        xmax <- max(ageDataOrig$depth)
        plot(ageData$depth,accRates,ylim=c(ymin,ymax),xlim=c(xmin,xmax),xlab="Depth [cm]",ylab="Rates [ka/cm]",cex.lab=4,cex.axis=4)
        if (Lsegments > 0) {
            for (i in seq(1,length(segPosPlot))) {
                par(new=TRUE)
                segPosDepth <- (ageData$depth[segPosPlot[i]]+ageData$depth[segPosPlot[i]+1])/2
                plot(c(segPosDepth,segPosDepth),c(ymin,ymax),ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=rgb(200/255,200/255,200/255),type="l",lwd="12",lty=2,xlab="",ylab="",xaxt="n",yaxt="n",cex=8)
            }
        }
        par(new=TRUE)
        polygon(c(ageData$depth, rev(ageData$depth)), c(fittedRates-fittedError,
                                                        rev(fittedRates+fittedError)), col=rgb(0.8,0.8,0.8),border=NA)
        par(new=TRUE)
        plot(ageData$depth,accRates,ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=rgb(42/255,190/255,230/255),pch=19,xlab="",ylab="",xaxt="n",yaxt="n",cex=8)
        par(new=TRUE)
        plot(ageData$depth,fittedRates,ylim=c(ymin,ymax),xlim=c(xmin,xmax),xlab="",pch=19,ylab="",xaxt="n",yaxt="n",col="red",cex=4)
        par(new=TRUE)
        plot(ageData$depth,fittedRates,ylim=c(ymin,ymax),xlim=c(xmin,xmax),lwd=8,xlab="",ylab="",xaxt="n",yaxt="n",type="l",col="red",cex=8)

        arrows(ageData$depth, accRates-sigmaRates, ageData$depth, 
               accRates+sigmaRates, length=0.2, angle=90, code=3,col="black",lwd=6,cex=8)

        legend(xmax, ymax, xjust=1,legend=c("data","fit","sequences"),
               col=c(rgb(42/255,190/255,230/255),"red",rgb(200/255,200/255,200/255)), lwd=c(NaN,8,12), lty=c(1,1,2), pch=c(19,19,NaN), cex=c(6), pt.cex=c(8,4,8))

        dev.off()
    }
    ##--------------------------------------------------------------------------------------------##




    ##--------------------------------------------------------------------------------------------##
    ## plotting ages
    for (i in plot_format){
    if (i == "png"){
            png(filename=paste(path,token,"_ages.png",sep=""),width=2000,height=1400,pointsize=10)#,width=11,height=8)
        }
    if (i == "pdf"){
            pdf(file=paste(path,token,"_ages.pdf",sep=""),width=40,height=28)#,width=11,height=8)
        }


        par(mar=c(14,14,4,2),mgp=c(10,4,0),bg=rgb(248/255,250/255,252/255))


        ymin <- min(c(ageData$age-ageData$error,ageDataOrig$age),na.rm=TRUE)
        ymax <- max(c(ageData$age+ageData$error,ageDataOrig$age,ageFit+ageFitError),na.rm=TRUE)
        xmin <- 0 ##min(ageDataOrig$depth)
        xmax <- max(ageDataOrig$depth)
        plot(ageData$depth,ageData$age,ylim=c(ymin,ymax),xlim=c(xmin,xmax),xlab="Depth [cm]",ylab="Age [ka]",cex.lab=4,cex.axis=4)
        if (Lsegments > 0) {
            for (i in seq(1,length(segPosPlot))) {
                par(new=TRUE)
                segPosDepth <- (ageData$depth[segPosPlot[i]]+ageData$depth[segPosPlot[i]+1])/2
                plot(c(segPosDepth,segPosDepth),c(ymin,ymax),ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=rgb(200/255,200/255,200/255),type="l",lwd="12",lty=2,xlab="",ylab="",xaxt="n",yaxt="n",cex=8)
            }
        }
        par(new=TRUE)
        polygon(c(ageData$depth, rev(ageData$depth)), c(ageFit-ageFitError,
                                                        rev(ageFit+ageFitError)), col=rgb(0.8,0.8,0.8),border=NA,xlim=c(xmin,xmax),)
        par(new=TRUE)
        plot(ageData$depth,ageData$age,ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=rgb(42/255,190/255,230/255),pch=19,xlab="",ylab="",xaxt="n",yaxt="n",cex=8)
        par(new=TRUE)
        plot(ageDataOrig$depth[negRatesIndex],ageDataOrig$age[negRatesIndex],ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=rgb(230/255,190/255,42/255),pch=19,xlab="",ylab="",xaxt="n",yaxt="n",cex=8)
        par(new=TRUE)
        plot(ageData$depth,ageFit,ylim=c(ymin,ymax),xlim=c(xmin,xmax),xlab="",pch=19,ylab="",xaxt="n",yaxt="n",col="red",cex=4)
        par(new=TRUE)
        plot(ageData$depth,ageFit,ylim=c(ymin,ymax),xlim=c(xmin,xmax),lwd=8,xlab="",ylab="",xaxt="n",yaxt="n",type="l",col="red",cex=8)

        arrows(ageDataOrig$depth, ageDataOrig$age-ageDataOrig$error, ageDataOrig$depth, 
               ageDataOrig$age+ageDataOrig$error, length=0.2, angle=90, code=3,col="black",lwd=6,cex=8)

        legend(xmax, ymin, xjust=1,yjust=0,legend=c("data","outliers","fit","sequences"),
               col=c(rgb(42/255,190/255,230/255),rgb(230/255,190/255,42/255),"red",rgb(200/255,200/255,200/255)), lwd=c(NaN,NaN,8,12), lty=c(1,1,1,2), pch=c(19,19,19,NaN), cex=c(6), pt.cex=c(8,8,4,NaN))


        dev.off()
    }
    ##--------------------------------------------------------------------------------------------##


    ##--------------------------------------------------------------------------------------------##
    ## plotting rates
    for (i in plot_format){
    if (i == "png"){
        png(filename=paste(path,token,"_rates.png",sep=""),width=2000,height=1400,pointsize=10)#,width=11,height=8)
        }
    if (i == "pdf"){
            pdf(file=paste(path,token,"_rates.pdf",sep=""),width=40,height=28)#,width=11,height=8)
        }
        
        par(mar=c(14,14,4,2),mgp=c(10,4,0),bg=rgb(248/255,250/255,252/255))
        
        ymin <- min(c(accRates-sigmaRates),na.rm=TRUE)
        ymax <- max(c(accRates+sigmaRates),na.rm=TRUE)

        xmin <- 0 ##min(ageDataOrig$depth)
        xmax <- max(ageDataOrig$depth)
        plot(ageData$depth,accRates,ylim=c(ymin,ymax),xlim=c(xmin,xmax),xlab="Depth [cm]",ylab="Rates [ka/cm]",cex.lab=4,cex.axis=4)
        if (Lsegments > 0) {
            for (i in seq(1,length(segPosPlot))) {
                par(new=TRUE)
                segPosDepth <- (ageData$depth[segPosPlot[i]]+ageData$depth[segPosPlot[i]+1])/2
                plot(c(segPosDepth,segPosDepth),c(ymin,ymax),ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=rgb(200/255,200/255,200/255),type="l",lwd="12",lty=2,xlab="",ylab="",xaxt="n",yaxt="n",cex=8)
            }
        }
        par(new=TRUE)
        polygon(c(ageData$depth, rev(ageData$depth)), c(BayesFittedRates-BayesSigmaRates,
                                                        rev(BayesFittedRates+BayesSigmaRates)), col=rgb(0.8,0.8,0.8),border=NA)
        par(new=TRUE)
        plot(ageData$depth,accRates,ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=rgb(42/255,190/255,230/255),pch=19,xlab="",ylab="",xaxt="n",yaxt="n",cex=8)
        par(new=TRUE)
        plot(ageData$depth,BayesFittedRates,ylim=c(ymin,ymax),xlim=c(xmin,xmax),xlab="",pch=19,ylab="",xaxt="n",yaxt="n",col="red",cex=4)
        par(new=TRUE)
        plot(ageData$depth,BayesFittedRates,ylim=c(ymin,ymax),xlim=c(xmin,xmax),lwd=8,xlab="",ylab="",xaxt="n",yaxt="n",type="l",col="red",cex=8)

        arrows(ageData$depth, accRates-sigmaRates, ageData$depth, 
               accRates+sigmaRates, length=0.2, angle=90, code=3,col="black",lwd=6,cex=8)

        legend(xmax, ymax, xjust=1,legend=c("data","fit","sequences"),
               col=c(rgb(42/255,190/255,230/255),"red",rgb(200/255,200/255,200/255)), lwd=c(NaN,8,12), lty=c(1,1,2), pch=c(19,19,NaN), cex=c(6), pt.cex=c(8,4,8))


        dev.off()
    }
    ##--------------------------------------------------------------------------------------------##



    ##--------------------------------------------------------------------------------------------##
    ## interpolate ages
    depthInterp <- seq(ageData$depth[1],ageData$depth[lengthData])
    lengthInterp <- length(depthInterp)
    ageFitInterp <- approx(x=ageData$depth,y=ageFit,xout=depthInterp,method="linear")
    ageFitErrorInterp <- approx(x=ageData$depth,y=ageFitError,xout=depthInterp,method="linear")

    ## interpolate rates
    ratesInterp <- approx(x=ageData$depth,y=BayesFittedRates,xout=depthInterp,method="linear")
    ratesErrorInterp <- approx(x=ageData$depth,y=BayesSigmaRates,xout=depthInterp,method="linear")


    ## write to disk
    out <- replicate(5,rep(0,lengthInterp))
    out[,1] <- ageFitInterp$x
    out[,2] <- ageFitInterp$y
    out[,3] <- ageFitErrorInterp$y
    out[,4] <- ratesInterp$y
    out[,5] <- ratesErrorInterp$y

    write("#depths ages sigma_ages rates sigma_rates",file=paste(path,token,"_interp.txt",sep=""),ncolumns=5)
    write(t(out),file=paste(path,token,"_interp.txt",sep=""),ncolumns=5,append=TRUE)

    

    ##fitted data
    out <- replicate(5,rep(0,lengthData))
    out[,1] <- ageData$depth
    out[,2] <- ageFit
    out[,3] <- ageFitError
    out[,4] <- BayesFittedRates
    out[,5] <- BayesSigmaRates

    
    write("#depths ages sigma_ages rates sigma_rates",file=paste(path,token,"_fit.txt",sep=""),ncolumns=5)
    write(t(out),file=paste(path,token,"_fit.txt",sep=""),ncolumns=5,append=TRUE)
    ##--------------------------------------------------------------------------------------------##

}
