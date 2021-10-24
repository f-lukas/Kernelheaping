#' @importFrom grDevices rainbow
#' @importFrom graphics box contour image legend lines par plot
#' @importFrom stats aggregate density dnorm optim pnorm qnorm quantile rgamma runif rnorm sd
#' @importFrom GB2 mlfit.gb2
#' @importFrom fitdistrplus fitdist
#' @importFrom sp point.in.polygon
#' @importFrom plyr round_any count
#' @importFrom magrittr "%>%"
#' @importFrom mvtnorm dmvnorm
#' @importFrom parallel clusterExport makeCluster stopCluster parLapply detectCores clusterEvalQ
#' @importFrom fastmatch fmatch

# create package environment for logging output
kernelheapingEnv <- new.env()
assign("logLevelInfo"  , TRUE ,  envir=kernelheapingEnv)
assign("logLevelDump"  , FALSE,  envir=kernelheapingEnv)
assign("logLevelTiming", FALSE,  envir=kernelheapingEnv)
assign("outputLevel"   , 1    ,  envir=kernelheapingEnv)

Rprx=function(rounds, Bias, RR, beta, gridx, ranefvalues, roundvalues){
  rprobs=rprobsRR2(RR,beta,gridx,ranefvalues)
  postMatrix=matrix(0,ncol=ncol(rprobs),nrow=2*nrow(rprobs))
  for(i in 1:length(rounds)){
    modulo=gridx%%rounds[i]
    lower=which(modulo < 0.5*rounds[i])
    higher=which(modulo > 0.5*rounds[i])
    postMatrix[i,lower]=Bias*rprobs[i,lower]
    postMatrix[i+length(rounds),higher]=(1-Bias)*rprobs[i,higher]
  }
  if(min(roundvalues)==0){postMatrix[c(1,1+length(rounds)),]=sum(postMatrix[c(1,1+length(rounds)),])/(2*length(gridx))}
  postMatrix=sweep(postMatrix,2,colSums(postMatrix),`/`, check.margin = FALSE) #\pi(R_i|X_i,\tau,a,\beta)
  return(postMatrix)
}

Rprx2=function(postMatrix, rounds, gridx){
  possW=seq(from=plyr::round_any(min(gridx),min(rounds)),to=round_any(max(gridx),min(rounds)),by=min(rounds))
  numbers=matrix(0,nrow=length(rounds)*2,ncol=length(possW),dimnames=list(NULL,possW))
  for(i in 1:length(c(rounds,rounds))){
    atemp=factor(plyr::round_any(gridx,accuracy=c(rounds,rounds)[i],f=round),levels=possW)
    numbers[i,]=tapply(postMatrix[i,],INDEX=list(atemp),FUN=function(x) sum(x,na.rm=T))
  }
  numbers[is.na(numbers)]=0
  sweep(numbers,2,colSums(numbers),`/`, check.margin = FALSE)# int \pi(R_i|W_i) dX_i
}

dbc <- function(gridx, x, bw, boundary=c(0, +Inf)){
  kernels <- sapply(x, function(y) {
    dens <- dnorm(gridx,mean = y,sd = bw)
    dens[gridx<boundary[1]|gridx>boundary[2]] <- 0
    dens
  })
  kernels <- sweep(kernels,2,colSums(kernels),`/`)
  rowMeans(kernels)/(gridx[2]-gridx[1])
}

rprobsRR2=function(RR,beta,gridx,ranefvalues=0){
  pgivenX=matrix(pnorm(rep(beta*log(abs(gridx))+ranefvalues,length(RR))+rep(RR,each=length(gridx))),
                 nrow=length(RR),byrow=T)
  diff(rbind(0,pgivenX,1))
}

logLik=function(par, new, rguesstemp, unequal, setBias, rounds,ranefvalues, roundvalues){
  RR=sort(par[1:(length(rounds)-1)])
  Bias=0
  beta=0
  if(setBias==TRUE) {Bias=par[length(rounds)]}
  if(unequal==TRUE) {beta=par[length(par)]}
  pgivenX=Rprx(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=beta, gridx=unique(new), ranefvalues=ranefvalues,
               roundvalues=roundvalues)+1E-96
  probs <- log(pgivenX[
    cumsum(rep(nrow(pgivenX),length.out=length(new)))-cumsum(((sequence(rle(new)$length)>1)*nrow(pgivenX)))-
      nrow(pgivenX)+rguesstemp])
  return(-sum(probs))
}

logLikRandomGamma=function(ranefvalues, RR, Bias, beta, new, rguesstemp, rounds, tau,roundvalues){
  pgivenX=Rprx(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=beta, gridx=unique(new), ranefvalues=ranefvalues,roundvalues=roundvalues)
  return(sum(log(pgivenX[cumsum(rep(nrow(pgivenX),length.out=length(new)))-
                           cumsum(((sequence(rle(new)$length)>1)*nrow(pgivenX)))-
                           nrow(pgivenX)+rguesstemp]))+sum(log(dnorm(ranefvalues,mean=0,sd=sqrt(tau)))))
}

#' Kernel density estimation for heaped data
#' @param xheaped heaped values from which to estimate density of x
#' @param rounds rounding values, numeric vector of length >=1
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param setBias if TRUE a rounding Bias parameter is estimated. For values above 0.5, the respondents
#' are more prone to round down, while for values < 0.5 they are more likely to round up
#' @param bw bandwidth selector method, defaults to "nrd0" see \code{density} for more options
#' @param boundary TRUE for positive only data (no positive density for negative values)
#' @param unequal if TRUE a probit model is fitted for the rounding probabilities with log(true value) as regressor
#' @param random if TRUE a random effect probit model is fitted for rounding probabilities
#' @param adjust as in \code{density}, the user can multiply the bandwidth by a certain factor such that bw=adjust*bw
#' @param weights optional numeric vector of sampling weights
#' @param recall if TRUE a recall error is introduced to the heaping model
#' @param recallParams recall error model parameters expression(nu) and expression(eta). Default is c(1/3, 1/3)
#' @return
#' The function returns a list object with the following objects (besides all input objects):
#' \item{\code{meanPostDensity}}{Vector of Mean Posterior Density}
#' \item{\code{gridx}}{Vector Grid on which density is evaluated}
#' \item{\code{resultDensity}}{Matrix with Estimated Density for each iteration}
#' \item{\code{resultRR}}{Matrix with rounding probability threshold values for each iteration (on probit scale)}
#' \item{\code{resultBias}}{Vector with estimated Bias parameter for each iteration}
#' \item{\code{resultBeta}}{Vector with estimated Beta parameter for each iteration}
#' \item{\code{resultX}}{Matrix of true latent values X estimates}
#' @examples
#' #Simple Rounding   ----------------------------------------------------------
#' xtrue=rnorm(3000)
#' xrounded=round(xtrue)
#' est <- dheaping(xrounded,rounds=1,burnin=20,samples=50)
#' plot(est,trueX=xtrue)
#'
#' #####################
#' #####Heaping
#' #####################
#' #Real Data Example  ----------------------------------------------------------
#' # Student learning hours per week
#' data(students)
#' xheaped <- as.numeric(na.omit(students$StudyHrs))
#' \dontrun{est <- dheaping(xheaped,rounds=c(1,2,5,10), boundary=TRUE, unequal=TRUE,burnin=20,samples=50)
#' plot(est)
#' summary(est)}
#'
#' #Simulate Data   ----------------------------------------------------------
#' Sim1 <- createSim.Kernelheaping(n=500, distribution="norm",rounds=c(1,10,100),
#' thresholds=c(-0.5244005, 0.5244005), sd=100)
#' \dontrun{est <- dheaping(Sim1$xheaped,rounds=Sim1$rounds)
#' plot(est,trueX=Sim1$x)}
#'
#' #Biased rounding
#' Sim2 <- createSim.Kernelheaping(n=500, distribution="gamma",rounds=c(1,2,5,10),
#'                      thresholds=c(-1.2815516, -0.6744898, 0.3853205),downbias=0.2,
#'                      shape=4,scale=8,offset=45)
#' \dontrun{est <- dheaping(Sim2$xheaped, rounds=Sim2$rounds, setBias=T, bw="SJ")
#' plot(est, trueX=Sim2$x)
#' summary(est)
#' tracePlots(est)}
#'
#' Sim3 <- createSim.Kernelheaping(n=500, distribution="gamma",rounds=c(1,2,5,10),
#' thresholds=c(1.84, 2.64, 3.05), downbias=0.75, Beta=-0.5, shape=4, scale=8)
#' \dontrun{est <- dheaping(Sim3$xheaped,rounds=Sim3$rounds,boundary=TRUE,unequal=TRUE,setBias=T)
#' plot(est,trueX=Sim3$x)}
#' @export
dheaping <- function(xheaped, rounds, burnin=5, samples=10, setBias=FALSE, weights=NULL,
                     bw= "nrd0", boundary=FALSE, unequal=FALSE, random=FALSE, adjust=1,
                     recall=F, recallParams=c(1/3,1/3)){
  if(is.null(weights)){weights=rep(1/length(xheaped),length(xheaped))}
  weights=weights/sum(weights)
  roundvalues=rounds
  if(min(rounds)==0){rounds[rounds==0]=0.5*min(rounds[rounds>0])}
  weights=weights[order(xheaped)]
  xheaped=xheaped[order(xheaped)] #round to lowest rounding value
  weights=weights[!is.na(xheaped)]
  xheaped=xheaped[!is.na(xheaped)]
  xheapedOriginal <- xheaped
  xheaped=plyr::round_any(xheaped,accuracy=min(rounds)) #round to lowest rounding value
  #Create grid
  stepsize=1/2^(which(diff(range(xheaped))/min(rounds)*2^(1:10)>100)[1]) #at least 100 grid points
  gridx=seq(min(xheaped)-max(rounds)*0.5+stepsize/2*min(rounds),
            max(xheaped)+max(rounds)*0.5-stepsize/2*min(rounds),
            by=min(rounds)*stepsize)
  if(boundary==TRUE|unequal==TRUE){gridx=gridx[gridx>0]}
  #Pilot Estimation
  if(boundary==FALSE){Mestimates <- density(xheaped,from=min(gridx),to=max(gridx),n=length(gridx),bw=max(rounds)*2,
                                            weights=weights)$y}
  if(boundary==TRUE){Mestimates <- dbc(gridx=gridx,x=xheaped,bw=max(rounds)*2)}
  #Starting values for x(new), Rounding values(rguesstemp), und Rounding thresholds(RR), Bias (Bias) and beta
  #rguesstemp<-sapply(1:length(xheaped),function(x) rounds[max(which(xheaped[x]%%rounds==0))])
  rguesstemp=rep(0,length(xheaped))
  heapedvalues=unique(xheaped)
  xtable<-table(xheaped)
  #Starting Values
  new<-xheaped
  RR<-qnorm(1:(length(rounds)-1)/(length(rounds)))
  Bias<-0
  beta<-0
  tau<-0.5
  ranefnames=unique(plyr::round_any(gridx,accuracy=min(rounds)))
  ranef=rnorm(length(ranefnames),sd=sqrt(tau))
  ranefvalues=rep(0,length(gridx))
  
  #Result matrices
  resultDensity=matrix(nrow=samples+burnin,ncol=length(gridx))
  resultX=matrix(nrow=samples+burnin,ncol=length(xheaped))
  resultRR=matrix(nrow=samples+burnin,ncol=length(rounds)-1)
  resultBias=matrix(nrow=samples+burnin,ncol=1)
  resultBeta=matrix(nrow=samples+burnin,ncol=1)
  resultRandom=matrix(nrow=samples+burnin,ncol=length(ranefnames),dimnames=list(1:(samples+burnin),ranefnames))
  resultTau=rep(0,times=samples+burnin)
  
  selectionGrid<-lapply(heapedvalues,function(k){
    selectionupper=lapply(1:length(rounds), function(x)
      which(gridx>k-rounds[x]*0.5&gridx<=k))
    selectionlower=lapply(1:length(rounds), function(x)
      which(gridx<k+rounds[x]*0.5&gridx>=k))
    selection=c(selectionlower,selectionupper)
  })
  
  for(j in 1:(burnin+samples)){
    Rprs1=Rprx(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=beta,
               gridx=gridx, ranefvalues=ranefvalues,roundvalues=roundvalues)
    Rprs2=Rprx2(Rprs1, rounds=rounds, gridx=gridx)[,as.character(heapedvalues)]
    if(recall==FALSE){
      getNew <- getNewX(heapedvalues, selectionGrid, Mestimates, Rprs1, Rprs2,
                        xtable, xheaped, gridx, xheapedOriginal, roundvalues, rguesstemp, new)
    }
    if(recall==TRUE){
      getNew <- getNewV(heapedvalues, selectionGrid, Mestimates, Rprs1, Rprs2, xtable, xheaped,
                        gridx, xheapedOriginal, roundvalues, rguesstemp, new, Bias, beta, ranefvalues, RR, recallParams)
    }
    new=getNew[[1]]
    rguesstemp=getNew[[2]]
    #ordered
    newx=getNew[[1]][order(getNew[[1]])]
    rguesstempx=getNew[[2]][order(getNew[[1]])]
    #Laplace-Approximation:
    par=RR
    if(setBias==TRUE){par=c(par,Bias)}
    if(unequal==TRUE){par=c(par,beta)}
    if(length(rounds)>1){
      ranefvaluestemp=ranefvalues[match(unique(newx),gridx)]
      laplace=optim(par,logLik,new=newx,rguesstemp=rguesstempx,
                    unequal=unequal, setBias=setBias,
                    ranefvalues=ranefvaluestemp,rounds=rounds ,hessian=T, roundvalues=roundvalues,
                    method="BFGS")
      par=MASS::mvrnorm(1,mu=laplace$par,Sigma=solve(laplace$hessian+diag(1E-7,length(par))))
      
      #     Metropolis-Hastings als Alternative
      #       parNew=par+rnorm(length(par),sd=0.05)
      #       alpha=min(1,exp(-logLik(par = parNew, new=new,rguesstemp=rguesstemp, unequal=unequal, setBias=setBias,
      #                              postMatrix=matrix(0,ncol=length(unique(new)),nrow=2*length(rounds),dimnames=list(NULL,unique(new[order(new)]))),
      #                              ranefvalues=ranefvaluestemp,rounds=rounds)+
      #                         logLik(par = par, new=new,rguesstemp=rguesstemp, unequal=unequal, setBias=setBias,
      #                                postMatrix=matrix(0,ncol=length(unique(new)),nrow=2*length(rounds),dimnames=list(NULL,unique(new[order(new)]))),
      #                                ranefvalues=ranefvaluestemp,rounds=rounds)
      #       ))
      #       if(runif(1)<alpha){par=parNew}
      
      if(setBias==TRUE) {Bias=par[length(rounds)]}
      if(unequal==TRUE) {beta=par[length(par)]}
      RR=sort(par[1:(length(rounds)-1)])
      if(random==TRUE){
        #Metropolis Hastings for RE
        ranefnamestemp=unique(plyr::round_any(newx,accuracy=min(rounds)))
        raneftemp=ranef[match(ranefnamestemp,ranefnames)]
        ranefnew=raneftemp
        for(k in 1:length(raneftemp)){
          ranefnew[k]=raneftemp[k]+rnorm(1,sd=0.2)
          alpha=min(1,exp(logLikRandomGamma(ranefvalues=ranefnew[k], RR=RR, Bias=Bias, beta=beta,
                                            new=newx[plyr::round_any(newx,accuracy=min(rounds))==ranefnamestemp[k]],
                                            rguesstemp=rguesstempx[plyr::round_any(newx,accuracy=min(rounds))==ranefnamestemp[k]],
                                            rounds=rounds, tau=tau, roundvalues=roundvalues)-
                            logLikRandomGamma(ranefvalues=raneftemp[k], RR=RR, Bias=Bias, beta=beta,
                                              new=newx[plyr::round_any(newx,accuracy=min(rounds))==ranefnamestemp[k]],
                                              rguesstemp=rguesstempx[plyr::round_any(newx,accuracy=min(rounds))==ranefnamestemp[k]],
                                              rounds=rounds, tau=tau ,roundvalues=roundvalues)
          ))
          if(is.na(alpha)){alpha=0}
          if(runif(1)<alpha){raneftemp[k]=ranefnew[k]}
        }
        tau=1/rgamma(1,shape=0.001+length(raneftemp)/2,scale=(0.001+0.5*sum((raneftemp)^2))^-1)
        ranef[match(ranefnamestemp,ranefnames)]=raneftemp
        ranefvalues=ranef[match(plyr::round_any(gridx,accuracy=min(rounds)),ranefnames)]
      }
    }
    if(0 %in% roundvalues){
      new[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]=xheapedOriginal[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]
    }
    if(recall==T){
      sdV <- recallParams[1]*c(rounds,rounds)[rguesstemp]+
        recallParams[2]*sapply(new,function(x)
          diff(c(0,pnorm(rep(beta*log(abs(x)),length(RR))+RR),1))%*%
            roundvalues)
      
      new=sapply(1:length(new), function(y) sample(x=gridx,size=1,
                                                   prob=dnorm(gridx,mean=new[y]+qnorm(pnorm(Bias))*sdV[y],
                                                              sd=sdV[y])*Mestimates))
    }
    
    h <- density(new,bw=bw)$bw*adjust
    if(boundary==TRUE){Mestimates <- dbc(gridx=gridx,x=new,bw=h)}
    if(boundary==FALSE){Mestimates <- density(new,from=min(gridx),to=max(gridx),n=length(gridx),bw=h,weights=weights)$y}
    #Save Results
    resultDensity[j,]=Mestimates
    resultRR[j,]=RR
    resultBias[j,]=pnorm(Bias)
    resultBeta[j,]=beta
    resultX[j,]=new
    resultTau[j]=tau
    resultRandom[j,]=ranef
    print(paste("Iteration:",j,"of", burnin+samples))
  }
  meanPostDensity=apply(resultDensity[-c(1:burnin),],2,mean)
  est<-list(meanPostDensity=meanPostDensity,resultDensity=resultDensity,resultRR=resultRR, resultBias=resultBias,
            resultBeta=resultBeta, resultX=resultX, xheaped=xheaped, gridx=gridx, boundary=boundary, rounds=roundvalues,
            setBias=setBias, unequal=unequal, bw=bw, burnin=burnin, samples=samples, adjust=adjust,
            resultTau=resultTau, resultRandom=resultRandom, random=random, weights=weights)
  class(est) <- "Kernelheaping"
  return(est)
}

#' Plot Kernel density estimate of heaped data naively and corrected by partly bayesian model
#' @param x Kernelheaping object produced by \code{dheaping} function
#' @param trueX optional, if true values X are known (in simulations, for example) the 'Oracle' density estimate is added as well
#' @param ... additional arguments given to standard plot function
#' @return plot with Kernel density estimates (Naive, Corrected and True (if provided))
#' @method plot Kernelheaping
#' @export
plot.Kernelheaping <- function(x,trueX=NULL, ...){
  if(x$boundary==FALSE){
    plot(density(x$xheaped,bw=x$bw,adjust=x$adjust, weights=x$weights),xlab="x",ylab="Density",main="", ...)
    if(!is.null(trueX)){lines(density(trueX,bw=x$bw,adjust=x$adjust, weights=x$weights),col="blue",lty=2,lwd=2)}
  }
  if(x$boundary==TRUE){
    plot(dbc(gridx=x$gridx,x=x$xheaped,bw=density(x$xheaped,bw=x$bw,adjust=x$adjust)$bw)~x$gridx,type="l",xlab="x",ylab="Density", main="", ...)
    if(!is.null(trueX)){lines(dbc(x$gridx, trueX, bw=density(trueX,bw=x$bw,adjust=x$adjust)$bw)~x$gridx,col="blue",lty=2,lwd=2)}
  }
  lines(x$meanPostDensity~x$gridx,col="red",lty=4,lwd=2)
  if(is.null(trueX)){legend("topright",c("Naive","Corrected"),col=c("black","red"),lty=c(1,4,2))}
  if(!is.null(trueX)){legend("topright",c("Naive","Corrected","Oracle"),col=c("black","red","blue"),lty=c(1,4,2),lwd=c(1,2,2))}
}

#' Prints some descriptive statistics (means and quantiles) for the estimated rounding, bias and acceleration (beta) parameters
#' @param object Kernelheaping object produced by \code{dheaping} function
#' @param ... unused
#' @return Prints summary statistics
#' @method summary Kernelheaping
#' @export
summary.Kernelheaping <- function(object, ...){
  Threshold=c(-Inf,colMeans(object$resultRR[-c(1:object$burnin),]))
  names(Threshold)=object$rounds
  cat("Mean Posterior Threshold Values:", "\n")
  cat(Threshold, "\n")
  if(object$setBias==TRUE){
    Bias=c(mean(object$resultBias[-c(1:object$burnin)]),quantile(object$resultBias[-c(1:object$burnin)],c(0.025,0.975)))
    names(Bias)[1]="Mean"
    cat("\n", "Bias, Mean and 95% credible interval:", "\n")
    cat(Bias, "\n")
  }
  if(object$unequal==TRUE){
    unequal=c(mean(object$resultBeta[-c(1:object$burnin)]),quantile(object$resultBeta[-c(1:object$burnin)],c(0.025,0.975)))
    names(unequal)[1]="Mean"
    cat("\n", "Beta, Mean and 95% credible interval:", "\n")
    cat(unequal, "\n")}
}

#' Plots some trace plots for the rounding, bias and acceleration (beta) parameters
#' @param x Kernelheaping object produced by \code{dheaping} function
#' @param ... additional arguments given to standard plot function
#' @return Prints summary statistics
#' @export
tracePlots <- function(x, ...){
  par(mfrow=c(2,ceiling(ncol(x$resultRR)/2)))
  for(i in 1:ncol(x$resultRR)){
    plot(x$resultRR[,i],xlab="iteration", ylab="value", main=paste("Traceplot Threshold", i), type="l")
  }
  par(mfrow=c(1,1))
  if(x$setBias==TRUE){
    plot(x$resultBias,xlab="iteration", ylab="value", main="Traceplot Bias", type="l")
  }
  
  if(x$unequal==TRUE){
    plot(x$resultBeta,xlab="iteration", ylab="value", main="Traceplot Beta", type="l")
  }
  if(x$random==TRUE){
    plot(x$resultTau,xlab="iteration", ylab="value", main="Traceplot tau", type="l")
  }
}

#' Simulation of heaping correction method
#' @param simRuns number of simulations runs
#' @param n sample size
#' @param distribution name of the distribution where random sampling is available, e.g. "norm"
#' @param rounds rounding values, numeric vector of length >=1
#' @param thresholds rounding thresholds
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param downbias Bias parameter used in the simulation
#' @param setBias if TRUE a rounding Bias parameter is estimated. For values above 0.5, the respondents
#' are more prone to round down, while for values < 0.5 they are more likely to round up
#' @param bw bandwidth selector method, defaults to "nrd0" see \code{density} for more options
#' @param offset location shift parameter used simulation in simulation
#' @param boundary TRUE for positive only data (no positive density for negative values)
#' @param Beta Parameter of the probit model for rounding probabilities used in simulation
#' @param unequal if TRUE a probit model is fitted for the rounding probabilities with log(true value) as regressor
#' @param adjust as in \code{density}, the user can multiply the bandwidth by a certain factor such that bw=adjust*bw
#' @param ... additional attributes handed over to \code{createSim.Kernelheaping}
#' @return List of estimation results
#' @examples
#' \dontrun{Sims1 <- sim.Kernelheaping(simRuns=2, n=500, distribution="norm",
#' rounds=c(1,10,100), thresholds=c(0.3,0.4,0.3), sd=100)}
#' @export
sim.Kernelheaping <- function(simRuns, n, distribution, rounds, thresholds, downbias=0.5, setBias = FALSE, Beta=0, unequal = FALSE, burnin = 5, samples = 10,
                              bw = "nrd0", offset=0, boundary = FALSE, adjust = 1, ...){
  lapply(1:simRuns, function(x){
    print(x)
    Sim <- createSim.Kernelheaping(n, distribution, rounds, thresholds, downbias=downbias, Beta=Beta, offset=offset, ...)
    est <- dheaping(Sim$xheaped, rounds=Sim$rounds, burnin = burnin, samples = samples, setBias = setBias,
                    bw = bw, boundary = boundary, unequal = unequal, adjust = adjust)
    
    if(est$boundary==T){
      naiveD  <- dbc(est$gridx, est$xheaped, density(est$xheaped,from=min(est$gridx),to=max(est$gridx),n=length(est$gridx),bw=est$bw,adjust=est$adjust)$bw)
      oracleD <- dbc(est$gridx, Sim$x, density(Sim$x,from=min(est$gridx),to=max(est$gridx),n=length(est$gridx),bw=est$bw,adjust=est$adjust)$bw)
    }
    if(est$boundary==F){
      naiveD  <- density(est$xheaped,from=min(est$gridx),to=max(est$gridx),n=length(est$gridx),bw=est$bw,adjust=est$adjust)$y
      oracleD <- density(Sim$x,from=min(est$gridx),to=max(est$gridx),n=length(est$gridx),bw=est$bw,adjust=est$adjust)$y
    }
    trueD <- do.call(paste("d",distribution,sep=""), list(est$grid-offset, ...))
    return(list(est=est,Sim=Sim,naiveD=naiveD,oracleD=oracleD,trueD=trueD))
  })
}
#' Simulation Summary
#' @param sim Simulation object returned from sim.Kernelheaping
#' @param coverage probability for computing coverage intervals
#' @return list with summary statistics
#' @export
simSummary.Kernelheaping <- function(sim, coverage=0.9){
  RMISE <- sqrt(sapply(sim,function(x) cbind(sum(diff(x$est$grid)[1]*(x$est$meanPostDensity-x$trueD)^2),
                                             sum(diff(x$est$grid)[1]*(x$naiveD-x$trueD)^2),
                                             sum(diff(x$est$grid)[1]*(x$oracleD-x$trueD)^2))))
  RRmeans <- sapply(sim,function(x) colMeans(as.matrix(x$est$resultRR[-c(1:x$est$burnin),],nrow=x$est$samples)))
  RRlower <- sapply(sim,function(x) apply(x$est$resultRR[-c(1:x$est$burnin),],2,quantile,0.5*(1-coverage)))
  RRupper <- sapply(sim,function(x) apply(x$est$resultRR[-c(1:x$est$burnin),],2,quantile,1-(0.5*(1-coverage))))
  Biasmean <- sapply(sim,function(x) mean(x$est$resultBias[-c(1:x$est$burnin),]))
  Biaslower <- sapply(sim,function(x) quantile(x$est$resultBias[-c(1:x$est$burnin),],0.5*(1-coverage)))
  Biasupper <- sapply(sim,function(x) quantile(x$est$resultBias[-c(1:x$est$burnin),],1-(0.5*(1-coverage))))
  Betamean <- sapply(sim,function(x) mean(x$est$resultBeta[-c(1:x$est$burnin),]))
  Betalower <- sapply(sim,function(x) quantile(x$est$resultBeta[-c(1:x$est$burnin),],0.5*(1-coverage)))
  Betaupper <- sapply(sim,function(x) quantile(x$est$resultBeta[-c(1:x$est$burnin),],1-(0.5*(1-coverage))))
  #coverage
  coverageRR <- rowSums(sapply(1:length(sim), function(x) sim[[x]]$Sim$thresholds > RRlower[,x] & sim[[x]]$Sim$thresholds < RRupper[,x]))/length(sim)
  coverageBias <- sum(sapply(1:length(sim), function(x) sim[[x]]$Sim$downbias>Biaslower[x] & sim[[x]]$Sim$downbias<Biasupper[x]))/length(sim)
  coverageBeta <- sum(sapply(1:length(sim), function(x) sim[[x]]$Sim$Beta>Betalower[x] & sim[[x]]$Sim$Beta<Betaupper[x]))/length(sim)
  return(list(RMISE=RMISE,RRmeans=RRmeans,RRlower=RRlower,RRupper=RRupper,coverageRR=coverageRR,
              Biasmean=Biasmean,Biaslower=Biaslower,Biasupper=Biasupper,coverageBias=coverageBias,
              Betamean=Betamean,Betalower=Betalower,Betaupper=Betaupper,coverageBeta=coverageBeta))
}
#' Create heaped data for Simulation
#' @param n sample size
#' @param distribution name of the distribution where random sampling is available, e.g. "norm"
#' @param rounds rounding values
#' @param thresholds rounding thresholds (for Beta=0)
#' @param offset certain value added to all observed random samples
#' @param downbias bias parameter
#' @param Beta acceleration paramter
#' @param ... additional attributes handed over to "rdistribution" (i.e. rnorm, rgamma,..)
#' @return List of heaped values, true values and input parameters
#' @export
createSim.Kernelheaping <- function(n, distribution, rounds, thresholds, offset=0, downbias=0.5, Beta=0, ...){
  xtrue=do.call(paste("r",distribution,sep=""), list(n, ...))+offset
  RR0=thresholds
  down=lapply(xtrue,function(x) (x%%rounds<=rounds*0.5)*downbias*rprobsRR2(RR=RR0,beta=Beta,gridx=x)+
                (x%%rounds>=rounds*0.5)*(1-downbias)*rprobsRR2(RR=RR0,beta=Beta,gridx=x))
  down=lapply(down, function(x) x/sum(x))
  roundings=sapply(down,function(z) sample(rounds,size=1,replace=TRUE,prob=z))
  xheaped=plyr::round_any(xtrue,accuracy=roundings, f=round)
  return(list(xheaped=xheaped, rounds=rounds, thresholds=thresholds, downbias=downbias, Beta=Beta, x=xtrue))
}
#' Bivariate kernel density estimation for rounded data
#' @param xrounded rounded values from which to estimate bivariate density, matrix with 2 columns (x,y)
#' @param roundvalue rounding value (side length of square in that the true value lies around the rounded one)
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param adaptive set to TRUE for adaptive bandwidth
#' @param gridsize number of evaluation grid points
#' @return
#' The function returns a list object with the following objects (besides all input objects):
#' \item{\code{Mestimates}}{kde object containing the corrected density estimate}
#' \item{\code{gridx}}{Vector Grid on which density is evaluated (x)}
#' \item{\code{gridy}}{Vector Grid on which density is evaluated (y)}
#' \item{\code{resultDensity}}{Array with Estimated Density for each iteration}
#' \item{\code{resultX}}{Matrix of true latent values X estimates}
#' \item{\code{delaigle}}{Matrix of Delaigle estimator estimates}
#' @examples
#' # Create Mu and Sigma  -----------------------------------------------------------
#' mu1 <- c(0, 0)
#' mu2 <- c(5, 3)
#' mu3 <- c(-4, 1)
#' Sigma1 <- matrix(c(4, 3, 3, 4), 2, 2)
#' Sigma2 <- matrix(c(3, 0.5, 0.5, 1), 2, 2)
#' Sigma3 <- matrix(c(5, 4, 4, 6), 2, 2)
#' # Mixed Normal Distribution -------------------------------------------------------
#' mus <- rbind(mu1, mu2, mu3)
#' Sigmas <- rbind(Sigma1, Sigma2, Sigma3)
#' props <- c(1/3, 1/3, 1/3)
#' \dontrun{xtrue=rmvnorm.mixt(n=1000, mus=mus, Sigmas=Sigmas, props=props)
#' roundvalue=2
#' xrounded=plyr::round_any(xtrue,roundvalue)
#' est <- dbivr(xrounded,roundvalue=roundvalue,burnin=5,samples=10)
#'
#' #Plot corrected and Naive distribution
#' plot(est,trueX=xtrue)
#' #for comparison: plot true density
#'  dens=dmvnorm.mixt(x=expand.grid(est$Mestimates$eval.points[[1]],est$Mestimates$eval.points[[2]]),
#'   mus=mus, Sigmas=Sigmas, props=props)
#'  dens=matrix(dens,nrow=length(est$gridx),ncol=length(est$gridy))
#'  contour(dens,x=est$Mestimates$eval.points[[1]],y=est$Mestimates$eval.points[[2]],
#'     xlim=c(min(est$gridx),max(est$gridx)),ylim=c(min(est$gridy),max(est$gridy)),main="True Density")}
#'@export
dbivr <- function(xrounded, roundvalue, burnin=2, samples=5, adaptive=FALSE, gridsize=200){
  #Create grid - sparr package requires identical grid length on both axis:
  gridx=seq(min(xrounded[,1])-0.5*roundvalue,max(xrounded[,1])+0.5*roundvalue,length=gridsize)
  gridy=seq(min(xrounded[,2])-0.5*roundvalue,max(xrounded[,2])+0.5*roundvalue,length=gridsize)
  
  #Pilot Estimation
  Mestimates <- ks::kde(x=xrounded, H=diag(c(roundvalue,roundvalue))^2,gridsize=c(length(gridx),length(gridy)),
                        xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),max(gridy)))
  
  #Result matrices
  resultDensity=array(dim=c(burnin+samples,length(gridx),length(gridy)))
  resultX=array(dim=c(samples+burnin,nrow(xrounded),2))
  
  rvalues=unique(xrounded) #unique rounded values in data
  selectionGrid<-lapply(1:nrow(rvalues),function(k){
    selectionX=which(Mestimates$eval.points[[1]]>=rvalues[k,1]-roundvalue*0.5&Mestimates$eval.points[[1]]<rvalues[k,1]+roundvalue*0.5)
    selectionY=which(Mestimates$eval.points[[2]]>=rvalues[k,2]-roundvalue*0.5&Mestimates$eval.points[[2]]<rvalues[k,2]+roundvalue*0.5)
    list(selectionX,selectionY)
  })
  
  
  #Delaigle Estimator
  delaigle=aggregate(list(length=rep(1,nrow(xrounded))), data.frame(xrounded), length)
  delaigle[,3]=delaigle[,3]/sum(delaigle[,3])/roundvalue^2
  delaigleest=matrix(0,nrow=length(gridx),ncol=length(gridy))
  for(i in 1:nrow(delaigle)){
    x=which(delaigle[i,1]==rvalues[,1]&delaigle[i,2]==rvalues[,2])
    delaigleest[selectionGrid[[x]][[1]],selectionGrid[[x]][[2]]]=delaigle[i,3]
  }
  
  #SEM Estimator
  for(j in 1:(burnin+samples)){
    new=c()
    for(i in 1:nrow(rvalues)){
      probs=as.vector(Mestimates$estimate[selectionGrid[[i]][[1]],selectionGrid[[i]][[2]]])
      points=cbind(rep(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]],
                       times=length(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]])),
                   rep(Mestimates$eval.points[[2]][selectionGrid[[i]][[2]]],
                       each=length(Mestimates$eval.points[[1]][selectionGrid[[i]][[1]]])))
      npoints=length(which(xrounded[,1]==rvalues[i,1]&xrounded[,2]==rvalues[i,2]))
      new=rbind(new,points[sample(1:nrow(points),size=npoints,replace=T,prob=probs),])
    }
    
    #recompute H
    if(adaptive==FALSE){
      H <- ks::Hpi(x = new, binned=TRUE) * 2
    }
    if(adaptive==TRUE){
      H <- ks::Hpi(x = new, binned=TRUE)
      H <- sqrt(sqrt(H[1,1]*H[2,2]))
    }
    #recompute density
    if(adaptive==FALSE){
      Mestimates <- ks::kde(x=new, H=H,gridsize=c(length(gridx),length(gridy)),bgridsize=c(length(gridx),length(gridy)),
                            xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),max(gridy)),binned=TRUE)
    }
    if(adaptive==TRUE){
      counts <- plyr::count(new)
      MestimatesAd <- sparr::bivariate.density(data=counts[,c(1:2)],pilotH=H,res=length(gridx),xrange=range(gridx),
                                               yrange=range(gridy),adaptive=TRUE,comment=FALSE,counts=counts[,3])
      Mestimates$estimate=MestimatesAd$Zm
    }
    Mestimates$estimate[is.na(Mestimates$estimate)]=1E-96
    resultDensity[j,,]=Mestimates$estimate
    resultX[j,,]=new
    print(paste("Iteration:",j,"of", burnin+samples))
  }
  Mestimates$estimate=apply(resultDensity[-c(1:burnin),,],c(2,3),mean)
  est<-list(Mestimates=Mestimates,resultDensity=resultDensity,resultX=resultX,
            xrounded=xrounded, gridx=gridx, gridy=gridy, roundvalue=roundvalue,
            burnin=burnin, samples=samples, adaptive=adaptive, delaigle=delaigleest)
  class(est) <- "bivrounding"
  return(est)
}

#' Plot Kernel density estimate of heaped data naively and corrected by partly bayesian model
#' @param x bivrounding object produced by \code{dbivr} function
#' @param trueX optional, if true values X are known (in simulations, for example) the 'Oracle' density estimate is added as well
#' @param ... additional arguments given to standard plot function
#' @return plot with Kernel density estimates (Naive, Corrected and True (if provided))
#' @method plot bivrounding
#' @export
plot.bivrounding <- function(x, trueX=NULL, ...){
  par(mfrow=c(2,2))
  
  #select Levels
  if(x$adaptive==FALSE){
    Naive=kde(x$xrounded,gridsize=c(length(x$gridx),length(x$gridy)),
              xmin=c(min(x$gridx),min(x$gridy)),xmax=c(max(x$gridx),max(x$gridy)))
    contour(x$Mestimates$estimate,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
            col=c("white",rainbow(6,start=0.5,end=0.6)),main="Corrected")
    contour(Naive$estimate,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
            col=c("white",rainbow(6,start=0.5,end=0.6)),main="Naive")
    if(!is.null(trueX)){
      Oracle=kde(trueX,gridsize=c(length(x$gridx),length(x$gridy)),
                 xmin=c(min(x$gridx),min(x$gridy)),xmax=c(max(x$gridx),max(x$gridy)))
      contour(Oracle$estimate,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
              col=c("white",rainbow(6,start=0.5,end=0.6)),main="Oracle")
    }
  }
  
  if(x$adaptive==TRUE){
    Naive=sparr::bivariate.density(data=x$xrounded,pilotH=ks::hpi(x = x$xrounded,binned=TRUE),
                                   res=length(x$gridx),xrange=range(x$gridx),
                                   yrange=range(x$gridy),adaptive=TRUE,comment=FALSE)
    contour(x$Mestimates$estimate,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
            col=c("white",rainbow(6,start=0.5,end=0.6)),main="Corrected")
    contour(Naive$Zm,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
            col=c("white",rainbow(6,start=0.5,end=0.6)),main="Naive")
    if(!is.null(trueX)){
      Oracle=sparr::bivariate.density(data=trueX,pilotH=ks::hpi(x = trueX),
                                      res=length(x$gridx),xrange=range(x$gridx),
                                      yrange=range(x$gridy),adaptive=TRUE,comment=FALSE)
      contour(Oracle$Zm,x=x$Mestimates$eval.points[[1]],y=x$Mestimates$eval.points[[2]],
              col=c("white",rainbow(6,start=0.5,end=0.6)),main="Oracle")
    }
  }
  
  image(x$delaigle,oldstyle=TRUE,x=x$Mestimates$eval.points[[1]],xlab="",ylab="",
        y=x$Mestimates$eval.points[[2]],col=c("white",rainbow(6,start=0.5,end=0.6)),
        main="Delaigle")
  box()
  par(mfrow=c(1,1))
}

getNewX <- function(heapedvalues, selectionGrid, Mestimates, Rprs1, Rprs2,
                    xtable, xheaped, gridx, xheapedOriginal, roundvalues, rguesstemp, new){
  for(i in 1:length(heapedvalues)){
    selection=selectionGrid[[i]]
    selectionprobs=lapply(1:length(selection),function(x)
      Mestimates[selection[[x]]]*
        #f(X_i,R_i|.)
        Rprs1[x,selection[[x]]]/(sum(Rprs1[x,selection[[x]]]))*
        (Rprs2[x,i]))
    selectionprobs <- unlist(selectionprobs)
    selectionprobs[is.na(selectionprobs)] <- 1e-16
    temprounds=unlist(lapply(1:length(selection), function(x) rep(x,times=length(selection[[x]]))))
    temp=sample(1:length(selectionprobs),size=xtable[i],prob=selectionprobs,replace=TRUE)
    rguesstemp[xheaped==heapedvalues[i]]=temprounds[temp]
    new[xheaped==heapedvalues[i]]=gridx[unlist(selection)[temp]]
  }
  if(0 %in% roundvalues){
    rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0] = sample(c(1,length(roundvalues)+1),
                                                                           length(rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]),replace=T)
    rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])==0&rguesstemp %in% c(1,length(roundvalues)+1)] <-
      rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])==0&rguesstemp %in% c(1,length(roundvalues)+1)]+1
  }
  return(list(new,rguesstemp))
}

getNewV <- function(heapedvalues, selectionGrid, Mestimates, Rprs1, Rprs2, xtable, xheaped,
                    gridx, xheapedOriginal, roundvalues, rguesstemp, new, Bias, beta, ranefvalues, RR, recallParams,rounds){
  for(i in 1:length(heapedvalues)){
    selection=selectionGrid[[i]]
    selectionprobs=lapply(1:length(selection),function(x)
      #f(X_i,R_i|.)
      Rprs1[x,selection[[x]]]/(sum(Rprs1[x,selection[[x]]]))*
        (Rprs2[x,i]))
    selectionprobs <- unlist(selectionprobs)
    selectionprobs[is.na(selectionprobs)] <- 1e-16
    
    temprounds=unlist(lapply(1:length(selection), function(x) rep(x,times=length(selection[[x]]))))
    
    sdV <- recallParams[1]*c(roundvalues,roundvalues)[temprounds]+
      recallParams[2]*sapply(unlist(selection),function(x)
        diff(c(0,pnorm(rep(beta*log(abs(gridx[x])),length(RR))+RR),1))%*%roundvalues)
    
    
    selectionprobs2=sapply(new[xheaped==heapedvalues[i]],function(x)
      dnorm(gridx[unlist(selection)],mean=x-qnorm(pnorm(Bias))*sdV,
            sd=sdV))
    
    temp=sapply(1:ncol(selectionprobs2), function(j)
      sample(1:length(selectionprobs),size=1,prob=selectionprobs*selectionprobs2[,j],
             replace=TRUE))
    
    rguesstemp[xheaped==heapedvalues[i]]=temprounds[temp]
    new[xheaped==heapedvalues[i]]=gridx[unlist(selection)[temp]]
  }
  if(0 %in% roundvalues){
    rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]=sample(c(1,length(roundvalues)+1),
                                                                           length(rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])!=0]),replace=T)
    rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])==0&rguesstemp %in% c(1,length(roundvalues)+1)] <-
      rguesstemp[xheapedOriginal%%min(roundvalues[roundvalues>0])==0&rguesstemp %in% c(1,length(roundvalues)+1)]+1
  }
  return(list(new,rguesstemp))
}

#' Kernel density estimation for classified data
#' @param xclass classified values; matrix with two columns: lower and upper value
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param boundary TRUE for positive only data (no positive density for negative values)
#' @param bw bandwidth selector method, defaults to "nrd0" see \code{density} for more options
#' @param evalpoints number of evaluation grid points
#' @param adjust as in \code{density}, the user can multiply the bandwidth by a certain factor such that bw=adjust*bw
#' @param dFunc character optional density (with "d", "p" and "q" functions) function name for parametric estimation such as "norm" "gamma" or "lnorm"
#' @return
#' The function returns a list object with the following objects (besides all input objects):
#' \item{\code{Mestimates}}{kde object containing the corrected density estimate}
#' \item{\code{gridx}}{Vector Grid on which density is evaluated}
#' \item{\code{resultDensity}}{Matrix with Estimated Density for each iteration}
#' \item{\code{resultX}}{Matrix of true latent values X estimates}
#' @examples
#' x=rlnorm(500, meanlog = 8, sdlog = 1)
#' classes <- c(0,500,1000,1500,2000,2500,3000,4000,5000,6000,8000,10000,15000,Inf)
#' xclass <- cut(x,breaks=classes)
#' xclass <- cbind(classes[as.numeric(xclass)], classes[as.numeric(xclass) + 1])
#' densityEst <- dclass(xclass=xclass, burnin=20, samples=50, evalpoints=1000)
#' plot(densityEst$Mestimates~densityEst$gridx ,lwd=2, type = "l")
#' @export
dclass <- function(xclass, burnin=2, samples=5, boundary=FALSE, bw="nrd0",
                   evalpoints=200, adjust=1, dFunc = NULL){
  if(max(xclass)==Inf){xclass[xclass == Inf]=3*max(xclass[xclass != Inf])}
  
  classmeans <- apply(xclass, 1, mean)
  
  #Create grid - sparr package requires identical grid length on both axis:
  gridx=seq(min(xclass),max(xclass),length=evalpoints)
  #Pilot Estimation
  if(is.null(dFunc)){
    if(boundary==FALSE){Mestimates <- density(classmeans,from=min(gridx),to=max(gridx),n=length(gridx),bw=2*max(xclass) / 10)$y}
    if(boundary==TRUE){Mestimates <- dbc(gridx=gridx,x=classmeans,bw=2*max(xclass) / 10)}
  } else {
    if(dFunc != "gb2"){
      pars <- fitdist(classmeans, dFunc, method = "mme")$estimate
    } else {
      pars <- mlfit.gb2(classmeans)[[2]]$par
    }
    args <- vector(mode = "list", length = 1 + length(pars))
    args[[1]] <- gridx
    args[2:(1+length(pars))] <- pars[1:length(pars)]
    Mestimates <- do.call(paste0("d", dFunc), args)
    Mestimates[is.na(Mestimates)] <- 0
    Mestimates[Mestimates == Inf] <- max(Mestimates[Mestimates != Inf])
  }
  
  #Result matrices
  resultDensity=matrix(ncol=c(burnin+samples),nrow=length(gridx))
  resultX=matrix(ncol=c(burnin+samples),nrow=length(xclass))
  
  if(!is.null(dFunc)){
    parsOut <- matrix(ncol=c(burnin+samples),nrow=length(pars))
  }
  
  selectionGrid<-lapply(1:nrow(xclass),function(k){
    selection=which(gridx>=xclass[k, 1]&gridx<xclass[k,2])
    selection})
  #SEM Estimator
  for(j in 1:(burnin+samples)){
    new=c()
    for(i in 1:nrow(xclass)){
      probs=as.vector(Mestimates[selectionGrid[[i]]])
      points=gridx[selectionGrid[[i]]]
      new=c(new,points[sample(1:length(points),size=1,replace=T,prob=probs)])
    }
    if(is.null(dFunc)){
      NewDensity <- density(new,from=min(gridx),to=max(gridx),n=length(gridx),bw=bw,adjust=adjust)
      Mestimates <- NewDensity$y
      if(boundary==TRUE){
        Mestimates <- dbc(gridx=gridx,x=new,bw=NewDensity$bw)
      }
    } else {
      if(dFunc != "gb2"){
        pars <- fitdist(new, dFunc, method = "mme")$estimate
      } else {
        pars <- mlfit.gb2(new)[[2]]$par
      }
      args <- vector(mode = "list", length = 1 + length(pars))
      args[[1]] <- gridx
      args[2:(1+length(pars))] <- pars[1:length(pars)]
      Mestimates <- do.call(paste0("d", dFunc), args)
      Mestimates[is.na(Mestimates)] <- 0
      
      Mestimates[Mestimates == Inf] <- max(Mestimates[Mestimates != Inf])
    }
    resultDensity[,j]=Mestimates
    resultX[,j]=new
    if(!is.null(dFunc)){
      parsOut[,j] <- pars
    }
    
    print(paste("Iteration:",j,"of", burnin+samples))
  }
  Mestimates=apply(resultDensity[,-c(1:burnin)],c(1),mean)
  if(is.null(dFunc)){
    est<-list(Mestimates=Mestimates,resultDensity=resultDensity,resultX=resultX,
              xclass=xclass, gridx=gridx,
              burnin=burnin, samples=samples)
  } else {
    est<-list(Mestimates=Mestimates,resultDensity=resultDensity,resultX=resultX,
              parsOut = parsOut,
              xclass=xclass, gridx=gridx,
              burnin=burnin, samples=samples)
    
  }
  class(est) <- "classDensity"
  return(est)
}

#' Bivariate Kernel density estimation for data classified in polygons or shapes
#' @param data data.frame with 3 columns: x-coordinate, y-coordinate (i.e. center of polygon) and number of observations in area.
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param adaptive TRUE for adaptive kernel density estimation
#' @param shapefile shapefile with number of polygons equal to nrow(data)
#' @param gridsize number of evaluation grid points
#' @param deleteShapes shapefile containing areas without observations
#' @param boundary boundary corrected kernel density estimate?
#' @param fastWeights if TRUE weigths for boundary estimation are only computed for first 10 percent of samples to speed up computation
#' @param numChains number of chains of SEM algorithm
#' @param numThreads number of threads to be used (only applicable if more than one chains)
#' @return
#' The function returns a list object with the following objects (besides all input objects):
#' \item{\code{Mestimates}}{kde object containing the corrected density estimate}
#' \item{\code{gridx}}{Vector Grid of x-coordinates on which density is evaluated}
#' \item{\code{gridy}}{Vector Grid of y-coordinates on which density is evaluated}
#' \item{\code{resultDensity}}{Matrix with Estimated Density for each iteration}
#' \item{\code{resultX}}{Matrix of true latent values X estimates}
#' @examples
#' \dontrun{
#' library(maptools)
#' 
#' # Read Shapefile of Berlin Urban Planning Areas (download available from:
#' # https://www.statistik-berlin-brandenburg.de/opendata/RBS_OD_LOR_2015_12.zip)
#' Berlin <- rgdal::readOGR("X:/SomeDir/RBS_OD_LOR_2015_12.shp") #(von daten.berlin.de)
#' 
#' # Get Dataset of Berlin Population (download available from:
#' # https://www.statistik-berlin-brandenburg.de/opendata/EWR201512E_Matrix.csv)
#' data <- read.csv2("X:/SomeDir/EWR201512E_Matrix.csv")
#' 
#' # Form Dataset for Estimation Process
#' dataIn <- cbind(t(sapply(1:length(Berlin@polygons),
#'  function(x) Berlin@polygons[[x]]@labpt)), data$E_E65U80)
#' 
#' #Estimate Bivariate Density
#' Est <- dshapebivr(data = dataIn, burnin = 5, samples = 10, adaptive = FALSE,
#'                  shapefile = Berlin, gridsize = 325, boundary = TRUE)}
#' 
#' # Plot Density over Area:
#' \dontrun{breaks <- seq(1E-16,max(Est$Mestimates$estimate),length.out = 20)
#' image.plot(x=Est$Mestimates$eval.points[[1]],y=Est$Mestimates$eval.points[[2]],
#'           z=Est$Mestimates$estimate, asp=1, breaks = breaks,
#'           col =  colorRampPalette(brewer.pal(9,"YlOrRd"))(length(breaks)-1))
#' plot(Berlin, add=TRUE)}
#' @export
dshapebivr <- function(data, burnin = 2, samples = 5, adaptive = FALSE, shapefile,
                       gridsize = 200, boundary = FALSE, deleteShapes = NULL,
                       fastWeights = TRUE, numChains = 1, numThreads=1){
  ###########################################################
  ##### get polygon shape coordinates from data #####
  ###########################################################
  
  npoints <- data[ ,3] #delete duplicated points
  
  if (length(gridsize) == 1) {  gridsize = c(gridsize, gridsize)  }
  
  #Create grid - sparr package requires identical grid length on both axis:
  gridx <- seq(shapefile@bbox[1,1],shapefile@bbox[1,2],length=gridsize[1])
  gridy <- seq(shapefile@bbox[2,1],shapefile@bbox[2,2],length=gridsize[2])
  grid <- as.matrix(expand.grid(gridx,gridy))
  
  #Pilot Estimation
  Mestimates <- ks::kde(x=data[,c(1,2)], H=diag(c(diff(shapefile@bbox[1,])/sqrt(length(shapefile@polygons)),
                                                  c(diff(shapefile@bbox[2,])/sqrt(length(shapefile@polygons)))))^2,
                        gridsize=c(length(gridx),length(gridy)),
                        xmin=c(min(gridx),min(gridy)),
                        xmax=c(max(gridx),max(gridy)), w = nrow(data)*data[,3]/sum(data[,3]))
  #Result matrices
  resultDensity=array(dim=c(numChains, burnin+samples,length(gridx),length(gridy)))
  resultX=array(dim=c(numChains, samples+burnin,sum(npoints),2))
  
  selectionGrid <- lapply(1:nrow(data),function(i){
    lapply(1:length(shapefile@polygons[[i]]@Polygons), function(j){
      if(shapefile@polygons[[i]]@Polygons[[j]]@hole == FALSE){
        which(point.in.polygon(grid[,1],grid[,2],shapefile@polygons[[i]]@Polygons[[j]]@coords[,1],
                               shapefile@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
      }
    }) %>% unlist
  })
  
  selectionGridHole <- lapply(1:nrow(data),function(i){
    lapply(1:length(shapefile@polygons[[i]]@Polygons), function(j){
      if(shapefile@polygons[[i]]@Polygons[[j]]@hole == TRUE){
        which(point.in.polygon(grid[,1],grid[,2],
                               shapefile@polygons[[i]]@Polygons[[j]]@coords[,1],
                               shapefile@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
      }
    }) %>% unlist
  })
  selectionGrid <- lapply(1:length(selectionGrid), function(x)
    setdiff(selectionGrid[[x]],selectionGridHole[[x]]))
  
  outside <- grid[-unlist(selectionGrid), ]
  
  unselectionGrid = NULL
  unselectionGridHoles = NULL
  if(!is.null(deleteShapes)){
    unselectionGrid <- lapply(1:length(deleteShapes@polygons),function(i){
      lapply(1:length(deleteShapes@polygons[[i]]@Polygons), function(j){
        if(deleteShapes@polygons[[i]]@Polygons[[j]]@hole == FALSE){
          which(point.in.polygon(grid[,1],grid[,2],
                                 deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,1],
                                 deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
        }
      }
      ) %>% unlist
    }) %>% unlist
    unselectionGridHoles <- lapply(1:length(deleteShapes@polygons),function(i){
      lapply(1:length(deleteShapes@polygons[[i]]@Polygons), function(j){
        if(deleteShapes@polygons[[i]]@Polygons[[j]]@hole == TRUE){
          which(point.in.polygon(grid[,1],grid[,2],
                                 deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,1],
                                 deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
        }
      }
      ) %>% unlist
    }) %>% unlist
    unselectionGrid <- setdiff(unselectionGrid, unselectionGridHoles)
    outside <- grid[-setdiff(unlist(selectionGrid), unselectionGrid), ]
    selectionGrid <- lapply(selectionGrid, function(x) setdiff(x,unselectionGrid))
  }
  
  
  inside <- grid[unlist(selectionGrid), ]
  
  #if no match in grid for shape -> get nearest match (possible reason grid too small..)
  for(j in which(sapply(selectionGrid,length)==0)){
    selectionGrid[[j]] <-  which.min(colSums((t(grid) - as.numeric(data[j,c(1,2)]))^2))
  }
  
  # decide wether to use threads or not
  numCoresToUse = min(numThreads, numChains, detectCores())
  if (numCoresToUse == 1)
  {
    for(c in seq(1,numChains))
    {
      res <- dshapebivr_calcChain(c,NA,Mestimates,burnin,samples,grid,gridx,gridy,
                                  selectionGrid,shapefile,npoints,adaptive,
                                  boundary,fastWeights,data,inside,outside,
                                  deleteShapes,unselectionGrid)
      
      resultDensity[c,,,]   <- res$resultDensity
      resultX[c,,,]         <- res$resultX
      Mestimates            <- res$Mestimates
    }
  } else
  {
    cl <- makeCluster(numCoresToUse, type="PSOCK", outfile="")
    # export necessary functions and libraries
    clusterEvalQ(cl, { library(ks); library(mvtnorm); library(dplyr); library(fastmatch) })
    
    # start parallel execution
    baseFilename = paste(tempdir(), "/", runif(1, min=1000000, max=9999999), "_", sep="")
    parLapply(cl, 1:numChains, dshapebivr_calcChain, baseFilename,
              Mestimates,burnin,samples,grid,gridx,gridy,
              selectionGrid,shapefile,npoints,adaptive,boundary,
              fastWeights,data,inside,outside,deleteShapes,unselectionGrid)
    stopCluster(cl)
    gc()
    
    # copy results to output data structures
    for(c in seq(1,numChains))
    {
      ret <- NULL
      load(paste(baseFilename, c, sep = ""))
      
      resultDensity[c,,,]   <- ret$resultDensity
      if (isOutputlevelHigh()) { resultX[c,,,]         <- ret$resultX }
      Mestimates            <- ret$Mestimates
      
      # resultDensity[c,,,]   <- clusterResults[[c]]$resultDensity
      # resultX[c,,,]         <- clusterResults[[c]]$resultX
      # Mestimates            <- clusterResults[[c]]$Mestimates
      
      rm(ret); gc()
    }
  }
  
  
  
  
  
  indexGridColumns <- c( length( dim(resultDensity[,-c(1:burnin),,]) ) -1,
                         length( dim(resultDensity[,-c(1:burnin),,]) ))
  Mestimates$estimate=apply(resultDensity[,-c(1:burnin),,],indexGridColumns,mean)
  est<-list(Mestimates=Mestimates,resultDensity=resultDensity,
            data=data, gridx=gridx, gridy=gridy, shapefile=shapefile,
            burnin=burnin, samples=samples, adaptive=adaptive)
  if (isOutputlevelHigh()) { est[["resultX"]]=resultX }
  class(est) <- "bivshape"
  return(est)
}



dshapebivr_calcChain <- function(chain,
                                 saveAsBaseFileName,
                                 Mestimates,
                                 burnin,
                                 samples,
                                 grid,
                                 gridx,
                                 gridy,
                                 selectionGrid,
                                 shapefile,
                                 npoints,
                                 adaptive,
                                 boundary,
                                 fastWeights,
                                 data,
                                 inside,
                                 outside,
                                 deleteShapes,
                                 unselectionGrid)
{
  printInfo("Start Chain ", chain)
  
  # return object, if saveAsBaseFileName is != NA it will be saved as RData File instead
  ret = list()
  ret$resultDensity=array(dim=c(burnin+samples,length(gridx),length(gridy)))
  ret$resultX=array(dim=c(samples+burnin,sum(npoints),2))
  
  #pre calculations for faster match usage (fmatch)
  Max_2SingleID = max(grid)
  inside_2SingleID = apply(inside, 1, function(x) { x[1]*Max_2SingleID + x[2] } )
  
  #SEM Estimator
  for(j in 1:(burnin+samples)){
    new=matrix(nrow=sum(npoints), ncol=2)
    newCnt=1
    printTiming("Calc Probabilities", {
      for(i in 1:length(selectionGrid)){
        probs=Mestimates$estimate[cbind(fmatch(grid[selectionGrid[[i]],1],
                                               Mestimates$eval.points[[1]]),
                                        fmatch(grid[selectionGrid[[i]],2],
                                               Mestimates$eval.points[[2]]))]+1E-10
        points=matrix(ncol=2,grid[selectionGrid[[i]],])
        if(length(selectionGrid[[i]])==0){
          points <- matrix(ncol=2,shapefile@polygons[[i]]@labpt)
          probs <- 1
        }
        if(npoints[i] > 0){
          if(npoints[i] == 1){
            sampleProp <- 1
          } else {
            sampleProp = sample(1:nrow(points),size=max(0,npoints[i],na.rm=T),replace=T,prob=probs)
          }
          new[newCnt:(newCnt+npoints[i]-1), ] = points[sampleProp,]
          newCnt = newCnt+npoints[i]
        }
      }
    })
    #recompute H
    printTiming("Recompute H", {
      if(adaptive==FALSE){
        H <- ks::Hpi(x = new, binned=TRUE)
      }
      if(adaptive==TRUE){
        H <- ks::Hpi(x = new, binned=TRUE)
        H <- sqrt(sqrt(H[1,1]*H[2,2]))
      }
    })
    
    w = rep(1, nrow(new))
    
    if(boundary == TRUE){
      printTiming("Calc Weights", {
        if(fastWeights == FALSE || j <= ceiling(0.1*(burnin+samples))){
          weights <- calcWeights_fast( inside  = inside ,
                                       outside = outside,
                                       gridx   = gridx  ,
                                       gridy   = gridy  ,
                                       H       = H         )
          
        }
      })
      
      printTiming("Match Weights", {
        new_2SingleID    = apply(new, 1, function(x) { x[1]*Max_2SingleID + x[2] } )
        w <- weights[ fmatch(new_2SingleID, inside_2SingleID, nomatch = 1) ]
      })
    }
    #recompute density
    if(adaptive==FALSE){
      Mestimates <- ks::kde(x=new, H=H,gridsize=c(length(gridx),length(gridy)),
                            bgridsize=c(length(gridx),length(gridy)),
                            xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),
                                                                 max(gridy)),
                            binned=TRUE, w = w/mean(w))
    }
    if(adaptive==TRUE){
      counts <- plyr::count(new)
      MestimatesAd <- sparr::bivariate.density(data=counts[,c(1:2)],pilotH=H,
                                               res=length(gridx),xrange=range(gridx),
                                               yrange=range(gridy),adaptive=TRUE,
                                               comment=FALSE, counts=counts[,3])
      Mestimates$estimate=MestimatesAd$Zm
    }
    if(!is.null(deleteShapes)){
      Mestimates$estimate[-setdiff(unlist(selectionGrid), unselectionGrid)] = 0
    } else{
      Mestimates$estimate[-(unlist(selectionGrid))] = 0
    }
    Mestimates$estimate[is.na(Mestimates$estimate)]=1E-96
    ret$resultDensity[j,,] <- Mestimates$estimate
    ret$resultX[j,,] <- new
    printInfo("Iteration: ",j," of ", burnin+samples," in chain ", chain)
  }
  
  ret$Mestimates = Mestimates
  
  if (is.na(saveAsBaseFileName))
  {
    return( ret )
  }
  else
  {
    save(ret, file = paste(saveAsBaseFileName, chain, sep=""), envir = environment())
    return( NA )
  }
}




calcWeights_fast <- function(inside, outside, gridx, gridy, H) {
  # idea: Create grid/image twice as the size of normal grid and calculate
  #       densisty one time with the mean in the middle of the grid.
  #       Then use only indizes-calculation to get the elements which
  #       needs to be summarized.
  #       To get the 1D array position in the 2D image the calculation
  #       1D_index = y_index * imageSize_x + x_index
  #       is used (lowest index is imageSize_x+1, because R starts with index 1).
  
  sizeX  <- length(gridx)
  sizeY  <- length(gridy)
  sizeX2 <- sizeX * 2
  sizeY2 <- sizeY * 2
  
  # get rastersize
  rX <- gridx[2] - gridx[1]
  rY <- gridy[2] - gridy[1]
  
  # create grid with doubled size and same rastersize as orgiginal grid
  doubleGridInd <- as.matrix(expand.grid( seq(1,sizeX2), seq(1,sizeY2) ))
  doubleGrid    <- cbind( doubleGridInd[,1] * rX, doubleGridInd[,2] *rY   )
  
  # calculate densisty only one time on the double grid with mean at (sizeX * rX, sizeY * rY)
  pilot     <- dmvnorm(x = doubleGrid, mean = c(sizeX * rX, sizeY * rY), sigma = H)
  pilot_sum <- sum(pilot)
  
  # store results of density in an array which is accessible with 1D index
  pilot_Ind1D <- rep(0, sizeX2*(sizeY2+1))
  pilot_Ind1D[ (doubleGridInd[,2]) * sizeX2 + doubleGridInd[,1] ] = pilot
  
  # get 2D indizes for inside and outside grid points
  insideInd  <- cbind( fmatch( inside [,1] , gridx ), fmatch( inside [,2] , gridy ) )
  outsideInd <- cbind( fmatch( outside[,1] , gridx ),  fmatch( outside[,2] , gridy ) )

  # even faster but not readable
  term1 <- sizeX2 * ( insideInd[,2] + sizeY ) + sizeX + insideInd[,1]
  term2 <- insideInd[,1] + insideInd[,2] * sizeX2
  
  weights  <- sapply(1:nrow(inside), function(x){
    ind <- term1 - term2[x]
    pilot_sum / sum( pilot_Ind1D[ ind ])
  })
  
  return(weights)
}


#' Bivariate Kernel density estimation for data classified in polygons or shapes
#' @param data data.frame with 4 columns: x-coordinate, y-coordinate (i.e. center of polygon) and number of observations in area for partial population and number of observations for complete observations.
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param adaptive TRUE for adaptive kernel density estimation
#' @param shapefile shapefile with number of polygons equal to nrow(data)
#' @param gridsize number of evaluation grid points
#' @param boundary boundary corrected kernel density estimate?
#' @param fastWeights if TRUE weigths for boundary estimation are only computed for first 10 percent of samples to speed up computation
#' @param deleteShapes shapefile containing areas without observations
#' @param numChains number of chains of SEM algorithm
#' @param numThreads number of threads to be used (only applicable if more than one chains)
#' @examples
#' \dontrun{
#' library(maptools)
#' 
#' # Read Shapefile of Berlin Urban Planning Areas (download available from:
#'   https://www.statistik-berlin-brandenburg.de/opendata/RBS_OD_LOR_2015_12.zip)
#' Berlin <- rgdal::readOGR("X:/SomeDir/RBS_OD_LOR_2015_12.shp") #(von daten.berlin.de)
#' 
#' # Get Dataset of Berlin Population (download available from:
#' # https://www.statistik-berlin-brandenburg.de/opendata/EWR201512E_Matrix.csv)
#' data <- read.csv2("X:/SomeDir/EWR201512E_Matrix.csv")
#' 
#' # Form Dataset for Estimation Process
#' dataIn <- cbind(t(sapply(1:length(Berlin@polygons),
#' function(x) Berlin@polygons[[x]]@labpt)), data$E_E65U80, data$E_E)
#' 
#' #Estimate Bivariate Proportions (may take some minutes)
#' PropEst <- dshapebivrProp(data = dataIn, burnin = 5, samples = 20, adaptive = FALSE,
#' shapefile = Berlin, gridsize=325, numChains = 16, numThreads = 4)}
#' 
#' # Plot Proportions over Area:
#' \dontrun{
#' breaks <- seq(0,0.4,by=0.025)
#' image.plot(x=PropEst$Mestimates$eval.points[[1]],y=PropEst$Mestimates$eval.points[[2]],
#'           z=PropEst$proportion+1E-96, asp=1, breaks = breaks,
#'           col =  colorRampPalette(brewer.pal(9,"YlOrRd"))(length(breaks)-1))
#' plot(Berlin, add=TRUE)}
#' @export
dshapebivrProp <- function(data, burnin=2, samples=5, adaptive=FALSE,
                           shapefile, gridsize=200, boundary = FALSE,
                           deleteShapes = NULL, fastWeights = TRUE,
                           numChains=1, numThreads=1){
  ###########################################################
  ##### get polygon shape coordinates from data #####
  ###########################################################
  pol.x <- list()
  pol.y <- list()
  for (i in 1:length(shapefile@polygons)) {
    pol.x[[i]] <- shapefile@polygons[[i]]@Polygons[[1]]@coords[,1]
    pol.y[[i]] <- shapefile@polygons[[i]]@Polygons[[1]]@coords[,2]
  }
  
  npointsAll <- data[ ,4] #delete duplicated points
  npoints <- data[ ,3] #delete duplicated points
  
  if (length(gridsize) == 1) { gridsize = c(gridsize, gridsize) }
  #Create grid - sparr package requires identical grid length on both axis:
  gridx=seq(shapefile@bbox[1,1],shapefile@bbox[1,2],length=gridsize[1])
  gridy=seq(shapefile@bbox[2,1],shapefile@bbox[2,2],length=gridsize[2])
  grid <- as.matrix(expand.grid(gridx,gridy))
  
  #Pilot Estimation
  MestimatesAll <- ks::kde(x=data[,c(1,2)], H = 10 * diag(c(diff(shapefile@bbox[1,])/sqrt(length(shapefile@polygons)),
                                                     c(diff(shapefile@bbox[2,])/sqrt(length(shapefile@polygons)))))^2,
                           gridsize=c(length(gridx),length(gridy)),
                           xmin=c(min(gridx),min(gridy)),
                           xmax=c(max(gridx),max(gridy)),w=nrow(data)*data[,4]/sum(data[,4]))
  
  Mestimates <- ks::kde(x=data[,c(1,2)], H = 10 * diag(c(diff(shapefile@bbox[1,])/sqrt(length(shapefile@polygons)),
                                                  c(diff(shapefile@bbox[2,])/sqrt(length(shapefile@polygons)))))^2,
                        gridsize=c(length(gridx),length(gridy)),
                        xmin=c(min(gridx),min(gridy)),
                        xmax=c(max(gridx),max(gridy)),w=nrow(data)*data[,3]/sum(data[,3]))
  
  #Result matrices
  resultDensity=array(dim=c(numChains, burnin+samples,length(gridx),length(gridy)))
  resultX=array(dim=c(numChains,samples+burnin,sum(npoints),2))
  proportionArray=array(dim=c(numChains,burnin+samples,length(gridx),length(gridy)))
  
  printTiming("Calc selectionGrid", {
    selectionGrid <- lapply(1:nrow(data),function(i){
      lapply(1:length(shapefile@polygons[[i]]@Polygons), function(j){
        if(shapefile@polygons[[i]]@Polygons[[j]]@hole == FALSE){
          which(point.in.polygon(grid[,1],grid[,2],shapefile@polygons[[i]]@Polygons[[j]]@coords[,1],
                                 shapefile@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
        }
      }) %>% unlist
    })
  })
  
  printTiming("Calc selectionGridHole", {
    selectionGridHole <- lapply(1:nrow(data),function(i){
      lapply(1:length(shapefile@polygons[[i]]@Polygons), function(j){
        if(shapefile@polygons[[i]]@Polygons[[j]]@hole == TRUE){
          which(point.in.polygon(grid[,1],grid[,2],
                                 shapefile@polygons[[i]]@Polygons[[j]]@coords[,1],
                                 shapefile@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
        }
      }) %>% unlist
    })
  })
  
  printTiming("Recalc selectionGrid", {
    selectionGrid <- lapply(1:length(selectionGrid), function(x) setdiff(selectionGrid[[x]],selectionGridHole[[x]]))
    outside <- grid[-unlist(selectionGrid), ]
  })
  
  printTiming("Consider DeleteShapes", {
    unselectionGrid = NULL
    unselectionGridHoles = NULL
    if(!is.null(deleteShapes)){
      unselectionGrid <- lapply(1:length(deleteShapes@polygons),function(i){
        lapply(1:length(deleteShapes@polygons[[i]]@Polygons), function(j){
          if(deleteShapes@polygons[[i]]@Polygons[[j]]@hole == FALSE){
            which(point.in.polygon(grid[,1],grid[,2],
                                   deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,1],
                                   deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
          }
        }
        ) %>% unlist
      }) %>% unlist
      unselectionGridHoles <- lapply(1:length(deleteShapes@polygons),function(i){
        lapply(1:length(deleteShapes@polygons[[i]]@Polygons), function(j){
          if(deleteShapes@polygons[[i]]@Polygons[[j]]@hole == TRUE){
            which(point.in.polygon(grid[,1],grid[,2],
                                   deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,1],
                                   deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
          }
        }
        ) %>% unlist
      }) %>% unlist
      unselectionGrid <- setdiff(unselectionGrid, unselectionGridHoles)
      outside <- grid[-setdiff(unlist(selectionGrid), unselectionGrid), ]
      selectionGrid <- lapply(selectionGrid, function(x) setdiff(x,unselectionGrid))
    }
  })
  
  inside <- grid[unlist(selectionGrid), ]  #SEM Estimator
  
  #if no match in grid for shape -> get nearest match (possible reason grid too small..)
  for(j in which(sapply(selectionGrid,length)==0)){
    selectionGrid[[j]] <-  which.min(colSums((t(grid) - as.numeric(data[j,c(1,2)]))^2))
  }
  
  rm(selectionGridHole, pol.x, pol.y); gc()
  
  # decide wether to use threads or not
  numCoresToUse = min(numThreads, numChains, detectCores())
  if (numCoresToUse == 1)
  {
    for(c in seq(1,numChains))
    {
      res <- dshapebivrProp_calcChain(c,NA,MestimatesAll,Mestimates,burnin,
                                      samples,grid,gridx,gridy,selectionGrid,
                                      shapefile,npointsAll,
                                      npoints,adaptive,boundary,fastWeights,
                                      data,inside,outside,deleteShapes,
                                      unselectionGrid)
      
      resultDensity[c,,,]   <- res$resultDensity
      if (isOutputlevelHigh()) { resultX[c,,,]         <- res$resultX }
      proportionArray[c,,,] <- res$proportionArray
      Mestimates            <- res$Mestimates
      
      rm(res); gc()
    }
  } else
  {
    cl <- makeCluster(numCoresToUse , type="PSOCK", outfile="")
    # export necessary functions and libraries
    clusterEvalQ(cl, { library(ks); library(mvtnorm); library(dplyr); library(fastmatch); library(Kernelheaping) })
    
    # generate random filename which is used for result files for each thread
    baseFilename = paste(tempdir(), "/", runif(1, min=1000000, max=9999999), "_", sep="")
    # start parallel execution
    parLapply(cl, 1:numChains, dshapebivrProp_calcChain, baseFilename,
              MestimatesAll, Mestimates, burnin, samples, grid,
              gridx, gridy, selectionGrid, shapefile, npointsAll,
              npoints, adaptive,boundary,fastWeights, data, inside,
              outside, deleteShapes, unselectionGrid)
    stopCluster(cl)
    rm(cl)
    gc()
    
    # copy results to output data structures
    for(c in seq(1,numChains))
    {
      ret <- NULL
      load(paste(baseFilename, c, sep = ""))
      
      resultDensity[c,,,]   <- ret$resultDensity
      if (isOutputlevelHigh()) { resultX[c,,,]         <- ret$resultX }
      proportionArray[c,,,] <- ret$proportionArray
      Mestimates            <- ret$Mestimates
      
      # resultDensity[c,,,]   <- clusterResults[[c]]$resultDensity
      # resultX[c,,,]         <- clusterResults[[c]]$resultX
      # proportionArray[c,,,] <- clusterResults[[c]]$proportionArray
      # Mestimates            <- clusterResults[[c]]$Mestimates
      
      rm(ret); gc()
    }
  }
  
  indexGridColumns <- c( length( dim(resultDensity[,-c(1:burnin),,]) ) -1,
                         length( dim(resultDensity[,-c(1:burnin),,]) ))
  Mestimates$estimate=apply(resultDensity[,-c(1:burnin),,],indexGridColumns,mean)
  est<-list(Mestimates=Mestimates,resultDensity=resultDensity,
            data=data, gridx=gridx, gridy=gridy, shapefile=shapefile,
            burnin=burnin, samples=samples, adaptive=adaptive, numChains=numChains,
            proportions = apply(proportionArray[,-c(1:burnin),,],indexGridColumns,mean))
  if (isOutputlevelHigh()) { est[["resultX"]]=resultX }
  class(est) <- "bivshape"
  return(est)
}


dshapebivrProp_calcChain <- function(chain,
                                     saveAsBaseFileName,
                                     MestimatesAll,
                                     Mestimates,
                                     burnin,
                                     samples,
                                     grid,
                                     gridx,
                                     gridy,
                                     selectionGrid,
                                     shapefile,
                                     npointsAll,
                                     npoints,
                                     adaptive,
                                     boundary,
                                     fastWeights,
                                     data,
                                     inside,
                                     outside,
                                     deleteShapes,
                                     unselectionGrid)
{
  printInfo("Start Chain ", chain)
  # return object; if saveAsBaseFileName is != NA it will be saved as RData File instead
  ret = list()
  ret$resultDensity=array(dim=c(burnin+samples,length(gridx),length(gridy)))
  ret$resultX=array(dim=c(samples+burnin,sum(npoints),2))
  ret$proportionArray=array(dim=c(burnin+samples,length(gridx),length(gridy)))
  
  #pre calculations for faster match usage (fmatch)
  Max_2SingleID = max(grid)
  inside_2SingleID = apply(inside, 1, function(x) { x[1]*Max_2SingleID + x[2] } )
  
  for(j in 1:(burnin+samples)){
    newAll=matrix(nrow=sum(npointsAll), ncol=2)
    newAllCnt=1
    new=matrix(nrow=sum(npoints), ncol=2)
    newCnt=1
    print(j)

    printTiming("Calc Probabilities", {
      for(i in 1:length(selectionGrid)){
        probsAll=MestimatesAll$estimate[cbind(match(grid[selectionGrid[[i]],1],
                                                    MestimatesAll$eval.points[[1]]),
                                              match(grid[selectionGrid[[i]],2],
                                                    MestimatesAll$eval.points[[2]]))]+1E-16
        probs=Mestimates$estimate[cbind(match(grid[selectionGrid[[i]],1],
                                              Mestimates$eval.points[[1]]),
                                        match(grid[selectionGrid[[i]],2],
                                              Mestimates$eval.points[[2]]))]+1E-16
        probsAll = probsAll/sum(probsAll)
        probs = probs/sum(probs)
        probsPROP=(probs/probsAll)*(sum(data[i,3])/sum(data[i,4]))
        points=matrix(ncol=2,grid[selectionGrid[[i]],])
        if(length(selectionGrid[[i]])==0){points <- matrix(ncol=2,shapefile@polygons[[i]]@labpt)
        probs=1
        probsAll=1}
        if(npointsAll[i]>0){
          sampleAll=sample(1:nrow(points),size=max(0,npointsAll[i],na.rm=T),replace=T,prob=probsAll)
          newAll[newAllCnt:(newAllCnt+npointsAll[i]-1), ] = points[sampleAll,]
          newAllCnt = newAllCnt+npointsAll[i]
        }
        if(npoints[i]>0){
          if (length(sampleAll) == 1) {
            sampleProp = sampleAll
          } else {
            sampleProp <- sample(sampleAll,
                                 size = max(0, npoints[i], na.rm = T), prob = probsPROP[sampleAll])
          }
          new[newCnt:(newCnt+npoints[i]-1), ] = points[sampleProp,]
          newCnt = newCnt+npoints[i]
        }
      }
    })
    #recompute H
    printTiming("Recompute H", {
      if(adaptive==FALSE){
        H <- ks::Hpi(x = new, binned=TRUE) * 2
      }
      if(adaptive==TRUE){
        H <- ks::Hpi(x = new, binned=TRUE)
        H <- sqrt(sqrt(H[1,1]*H[2,2]))
      }
    })
    
    wAll = rep(1, nrow(newAll))
    w = rep(1, nrow(new))
    
    if(boundary == TRUE){
      printTiming("Calc Weights", {
        if(fastWeights == FALSE || j <= ceiling(0.1*(burnin+samples))){
          weights <- calcWeights_fast( inside  = inside ,
                                       outside = outside,
                                       gridx   = gridx  ,
                                       gridy   = gridy  ,
                                       H       = H         )
          
        }
      })
      
      printTiming("Match Weights", {
        new_2SingleID    = apply(new, 1, function(x) { x[1]*Max_2SingleID + x[2] } )
        newAll_2SingleID = apply(newAll, 1, function(x) { x[1]*Max_2SingleID + x[2] } )
        w <- weights[ fmatch(new_2SingleID, inside_2SingleID, nomatch = 1) ]
        wAll <- weights[ fmatch(newAll_2SingleID, inside_2SingleID, nomatch = 1) ]
      })
    }
    #recompute density
    printTiming("Recompute Density", {
      if(adaptive==FALSE){
        MestimatesAll <- ks::kde(x=newAll, H=H,gridsize=c(length(gridx),length(gridy)),bgridsize=c(length(gridx),length(gridy)),
                                 xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),max(gridy)),binned=TRUE, w = wAll/mean(wAll))
        Mestimates <- ks::kde(x=new, H=H,gridsize=c(length(gridx),length(gridy)),bgridsize=c(length(gridx),length(gridy)),
                              xmin=c(min(gridx),min(gridy)),xmax=c(max(gridx),max(gridy)),binned=TRUE, w = w/mean(w))
      }
      if(adaptive==TRUE){
        counts <- plyr::count(new)
        MestimatesAd <- sparr::bivariate.density(data=counts[,c(1:2)],pilotH=H,res=length(gridx),xrange=range(gridx),
                                                 yrange=range(gridy),adaptive=TRUE,
                                                 comment=FALSE, counts=counts[,3])
        Mestimates$estimate=MestimatesAd$Zm
      }
      #weird behaviour of ks::kde returning negative densities (probably due to numeric problems)
      Mestimates$estimate[Mestimates$estimate < 0] <- 0
      MestimatesAll$estimate[MestimatesAll$estimate < 0] <- 0
      
    })
    
    printTiming("Delete Shapes", {
      if(!is.null(deleteShapes)){
        Mestimates$estimate[-setdiff(unlist(selectionGrid), unselectionGrid)] = 0
        MestimatesAll$estimate[-setdiff(unlist(selectionGrid), unselectionGrid)] = 0
      } else{
        Mestimates$estimate[-(unlist(selectionGrid))] = 0
        MestimatesAll$estimate[-(unlist(selectionGrid))] = 0
      }
    })
    
    printTiming("Match Density", {
      densityAll <- matrix(NA,ncol=ncol(Mestimates$estimate), nrow=nrow(Mestimates$estimate))
      densityPart <- matrix(NA,ncol=ncol(Mestimates$estimate), nrow=nrow(Mestimates$estimate))
      densityAll[cbind(match(grid[selectionGrid %>% unlist ,1],MestimatesAll$eval.points[[1]]),
                       match(grid[selectionGrid %>% unlist,2],MestimatesAll$eval.points[[2]]))] =
        MestimatesAll$estimate[cbind(match(grid[selectionGrid %>% unlist ,1],MestimatesAll$eval.points[[1]]),
                                     match(grid[selectionGrid %>% unlist,2],MestimatesAll$eval.points[[2]]))]
      densityPart[cbind(match(grid[selectionGrid %>% unlist ,1],MestimatesAll$eval.points[[1]]),
                        match(grid[selectionGrid %>% unlist,2],MestimatesAll$eval.points[[2]]))] =
        Mestimates$estimate[cbind(match(grid[selectionGrid %>% unlist ,1],MestimatesAll$eval.points[[1]]),
                                  match(grid[selectionGrid %>% unlist,2],MestimatesAll$eval.points[[2]]))]
    })
    
    printTiming("Assignments", {
      densityPart <- densityPart/sum(densityPart, na.rm=TRUE)
      densityAll <- densityAll/sum(densityAll, na.rm=TRUE)
      densityAll <- densityAll+1E-96
      proportion <- (densityPart/densityAll) * (sum(data[,3])/sum(data[,4]))
      proportion[proportion<0]=0
      proportion[proportion>1]=1
      
      Mestimates$estimate[is.na(Mestimates$estimate)]=1E-96
      ret$resultDensity[j,,] <- Mestimates$estimate
      ret$proportionArray[j,,] <- proportion
      ret$resultX[j,,] <- new
    })
    printInfo("Iteration: ",j," of ", burnin+samples," in chain ", chain)
  }
  
  ret$Mestimates = Mestimates
  
  if (is.na(saveAsBaseFileName))
  {
    return( ret )
  }
  else
  {
    save(ret, file = paste(saveAsBaseFileName, chain, sep=""), envir = environment())
    return( NA )
  }
}

#' 3d Kernel density estimation for data classified in polygons or shapes
#' @param data data.frame with 5 columns: x-coordinate, y-coordinate (i.e. center of polygon) and number of observations in area for partial population and number of observations for complete observations and third variable.
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param shapefile shapefile with number of polygons equal to nrow(data) / length(unique(data[,5]))
#' @param gridsize number of evaluation grid points
#' @param boundary boundary corrected kernel density estimate?
#' @param fastWeights if TRUE weigths for boundary estimation are only computed for first 10 percent of samples to speed up computation
#' @param deleteShapes shapefile containing areas without observations
#' @param numChains number of chains of SEM algorithm
#' @param numThreads number of threads to be used (only applicable if more than one chains)
#' @export
dshape3dProp <- function(data, burnin=2, samples=5,
                         shapefile, gridsize=200, boundary = FALSE,
                         deleteShapes = NULL, fastWeights = TRUE,
                         numChains=1, numThreads=1){
  ###########################################################
  ##### get polygon shape coordinates from data #####
  ###########################################################
  
  # adaptive not implemented
  adaptive <- FALSE
  
  # order data
  data <- data[order(data[,5]),]
  
  pol.x <- list()
  pol.y <- list()
  for (i in 1:length(shapefile@polygons)) {
    pol.x[[i]] <- shapefile@polygons[[i]]@Polygons[[1]]@coords[,1]
    pol.y[[i]] <- shapefile@polygons[[i]]@Polygons[[1]]@coords[,2]
  }
  
  npointsAll <- data[ ,c(4,5)] # number of observations for complete observations and third variable
  npoints <- data[ ,c(3,5)] # number of observations in area for partial population and third variable
  
  if (length(gridsize) == 1) { gridsize = c(gridsize, gridsize) }
  #Create grid - sparr package requires identical grid length on both axis:
  gridx=seq(shapefile@bbox[1,1],shapefile@bbox[1,2],length=gridsize[1])
  gridy=seq(shapefile@bbox[2,1],shapefile@bbox[2,2],length=gridsize[2])
  grid <- as.matrix(expand.grid(gridx,gridy))
  
  # calculate weights for for complete observations and each unique observation of third variable
  dataAllWeight <- sapply(unique(data[,5]), function(x){
    d <- data[data[,5]==x,]
    nrow(d)*d[,4]/sum(d[,4])
  })
  
  # calculate weights for for partial population and each unique observation of third variable
  dataWeight <- sapply(unique(data[,5]), function(x){
    d <- data[data[,5]==x,]
    nrow(d)*d[,3]/sum(d[,3])
  })
  
  #Pilot Estimation 
  # 3 dimensional Kernel Density Estimation
  # added new weight variable and adjusted x, H, gridsize, xmin and xmax for third dimension
  MestimatesAll <- ks::kde(x=data[,c(1,2,5)],H = 10 * diag(c(diff(shapefile@bbox[1,])/sqrt(length(shapefile@polygons)),
                                                             c(diff(shapefile@bbox[2,])/sqrt(length(shapefile@polygons))),
                                                             c(diff(c(min(data[,5]),max(data[,5]))))))^2,
                           gridsize=c(length(gridx),length(gridy),length(unique(data[,5]))),
                           xmin=c(min(gridx),min(gridy),min(data[,5])),
                           xmax=c(max(gridx),max(gridy),max(data[,5])),
                           w=unlist(as.vector(dataAllWeight)))
  
  Mestimates <- ks::kde(x=data[,c(1,2,5)],H = 10 * diag(c(diff(shapefile@bbox[1,])/sqrt(length(shapefile@polygons)),
                                                          c(diff(shapefile@bbox[2,])/sqrt(length(shapefile@polygons))),
                                                          c(diff(c(min(data[,5]),max(data[,5]))))))^2,
                        gridsize=c(length(gridx),length(gridy),length(unique(data[,5]))),
                        xmin=c(min(gridx),min(gridy),min(data[,5])),
                        xmax=c(max(gridx),max(gridy),max(data[,5])),
                        w=unlist(as.vector(dataWeight)))
  
  #Result matrices
  resultDensity=array(dim=c(numChains, burnin+samples,length(gridx),length(gridy),length(unique(npointsAll[,2]))))
  resultX=array(dim=c(numChains,samples+burnin,sum(npoints[,1]),3))
  proportionArray=array(dim=c(numChains, burnin+samples,length(gridx),length(gridy),length(unique(npointsAll[,2]))))
  
  printTiming("Calc selectionGrid", {
    selectionGrid <- lapply(1:nrow(unique(data[,1:2])),function(i){
      lapply(1:length(shapefile@polygons[[i]]@Polygons), function(j){
        if(shapefile@polygons[[i]]@Polygons[[j]]@hole == FALSE){
          which(point.in.polygon(grid[,1],grid[,2],shapefile@polygons[[i]]@Polygons[[j]]@coords[,1],
                                 shapefile@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
        }
      }) %>% unlist
    })
  })
  
  printTiming("Calc selectionGridHole", {
    selectionGridHole <- lapply(1:nrow(unique(data[,1:2])),function(i){
      lapply(1:length(shapefile@polygons[[i]]@Polygons), function(j){
        if(shapefile@polygons[[i]]@Polygons[[j]]@hole == TRUE){
          which(point.in.polygon(grid[,1],grid[,2],
                                 shapefile@polygons[[i]]@Polygons[[j]]@coords[,1],
                                 shapefile@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
        }
      }) %>% unlist
    })
  })
  
  printTiming("Recalc selectionGrid", {
    selectionGrid <- lapply(1:length(selectionGrid), function(x) setdiff(selectionGrid[[x]],selectionGridHole[[x]]))
    outside <- grid[-unlist(selectionGrid), ]
  })
  
  printTiming("Consider DeleteShapes", {
    unselectionGrid = NULL
    unselectionGridHoles = NULL
    if(!is.null(deleteShapes)){
      unselectionGrid <- lapply(1:length(deleteShapes@polygons),function(i){
        lapply(1:length(deleteShapes@polygons[[i]]@Polygons), function(j){
          if(deleteShapes@polygons[[i]]@Polygons[[j]]@hole == FALSE){
            which(point.in.polygon(grid[,1],grid[,2],
                                   deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,1],
                                   deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
          }
        }
        ) %>% unlist
      }) %>% unlist
      unselectionGridHoles <- lapply(1:length(deleteShapes@polygons),function(i){
        lapply(1:length(deleteShapes@polygons[[i]]@Polygons), function(j){
          if(deleteShapes@polygons[[i]]@Polygons[[j]]@hole == TRUE){
            which(point.in.polygon(grid[,1],grid[,2],
                                   deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,1],
                                   deleteShapes@polygons[[i]]@Polygons[[j]]@coords[,2])==1)
          }
        }
        ) %>% unlist
      }) %>% unlist
      unselectionGrid <- setdiff(unselectionGrid, unselectionGridHoles)
      outside <- grid[-setdiff(unlist(selectionGrid), unselectionGrid), ]
      selectionGrid <- lapply(selectionGrid, function(x) setdiff(x,unselectionGrid))
    }
  })
  
  inside <- grid[unlist(selectionGrid), ]  #SEM Estimator
  
  #if no match in grid for shape -> get nearest match (possible reason grid too small..)
  for(j in which(sapply(selectionGrid,length)==0)){
    selectionGrid[[j]] <-  which.min(colSums((t(grid) - as.numeric(data[j,c(1,2)]))^2))
  }
  
  rm(selectionGridHole, pol.x, pol.y); gc()
  
  # decide wether to use threads or not
  numCoresToUse = min(numThreads, numChains, detectCores())
  if (numCoresToUse == 1)
  {
    for(c in seq(1,numChains))
    {
      res <- dshape3dProp_calcChain(c,NA,MestimatesAll,Mestimates,burnin,
                                    samples,grid,gridx,gridy,selectionGrid,
                                    shapefile,npointsAll,
                                    npoints,adaptive,boundary,fastWeights,
                                    data,inside,outside,deleteShapes,
                                    unselectionGrid)
      
      resultDensity[c,,,,]   <- res$resultDensity
      if (isOutputlevelHigh()) { resultX[c,,,]         <- res$resultX }
      proportionArray[c,,,,] <- res$proportionArray
      Mestimates            <- res$Mestimates
      
      rm(res); gc()
    }
  } else
  {
    cl <- makeCluster(numCoresToUse , type="PSOCK", outfile="")
    # export necessary functions and libraries
    clusterEvalQ(cl, { library(ks); library(mvtnorm); library(dplyr); library(fastmatch); library(Kernelheaping) })
    
    # generate random filename which is used for result files for each thread
    baseFilename = paste(tempdir(), "/", runif(1, min=1000000, max=9999999), "_", sep="")
    # start parallel execution
    parLapply(cl, 1:numChains, dshape3dProp_calcChain, baseFilename,
              MestimatesAll, Mestimates, burnin, samples, grid,
              gridx, gridy, selectionGrid, shapefile, npointsAll,
              npoints, adaptive,boundary,fastWeights, data, inside,
              outside, deleteShapes, unselectionGrid)
    
    stopCluster(cl)
    rm(cl)
    gc()
    
    # copy results to output data structures
    for(c in seq(1,numChains))
    {
      ret <- NULL
      load(paste(baseFilename, c, sep = ""))
      
      resultDensity[c,,,,]   <- ret$resultDensity
      if (isOutputlevelHigh()) { resultX[c,,,]         <- ret$resultX }
      proportionArray[c,,,,] <- ret$proportionArray
      Mestimates            <- ret$Mestimates
      
      # resultDensity[c,,,]   <- clusterResults[[c]]$resultDensity
      # resultX[c,,,]         <- clusterResults[[c]]$resultX
      # proportionArray[c,,,] <- clusterResults[[c]]$proportionArray
      # Mestimates            <- clusterResults[[c]]$Mestimates
      
      rm(ret); gc()
    }
  }
  
  indexGridColumns <- c( length( dim(resultDensity[,-c(1:burnin),,,]) ) -2,
                         length( dim(resultDensity[,-c(1:burnin),,,]) ) -1,
                         length( dim(resultDensity[,-c(1:burnin),,,]) ))
  Mestimates$estimate=apply(resultDensity[,-c(1:burnin),,,],indexGridColumns,mean)
  est<-list(Mestimates=Mestimates,resultDensity=resultDensity,
            data=data, gridx=gridx, gridy=gridy, shapefile=shapefile,
            burnin=burnin, samples=samples, adaptive=adaptive, numChains=numChains,
            proportions = apply(proportionArray[,-c(1:burnin),,,],indexGridColumns,mean))
  if (isOutputlevelHigh()) { est[["resultX"]]=resultX }
  class(est) <- "3dshape"
  return(est)
}


dshape3dProp_calcChain <- function(chain,
                                   saveAsBaseFileName,
                                   MestimatesAll,
                                   Mestimates,
                                   burnin,
                                   samples,
                                   grid,
                                   gridx,
                                   gridy,
                                   selectionGrid,
                                   shapefile,
                                   npointsAll,
                                   npoints,
                                   adaptive,
                                   boundary,
                                   fastWeights,
                                   data,
                                   inside,
                                   outside,
                                   deleteShapes,
                                   unselectionGrid)
{
  printInfo("Start Chain ", chain)
  # return object; if saveAsBaseFileName is != NA it will be saved as RData File instead
  ret = list()
  ret$resultDensity=array(dim=c(burnin+samples,length(gridx),length(gridy),length(unique(npointsAll[,2]))))
  ret$resultX=array(dim=c(samples+burnin,sum(npoints[,1]),3))
  ret$proportionArray=array(dim=c(burnin+samples,length(gridx),length(gridy),length(unique(npointsAll[,2]))))
  
  #pre calculations for faster match usage (fmatch)
  Max_2SingleID = max(grid)
  inside_2SingleID = apply(inside, 1, function(x) { x[1]*Max_2SingleID + x[2] } )
  
  for(j in 1:(burnin+samples)){
    
    printInfo("Start Iteration: ",j," of ", burnin+samples," in chain ", chain)
    
    newAll=NULL
    new=NULL
    
    # added a new loop:
    # in dshapebivrProp() sampling is done seperately for each area (i) from the 2d density
    # in dshapebivrProp3d() sampling is done seperately for each time period (k) and each area (i) from the 3d density
    
    for(k in 1:length(unique(npointsAll[,2]))){
      
      print(paste0("... ",k," / ",length(unique(npointsAll[,2]))))
      
      # subset data and npoints, keeping only observation k of third variable
      # do steps similar as in 2d case and later bind together sampled geocoordinates for each time period 
      datak <- data[data[,5]==unique(data[,5])[k],]
      
      npointsAllk <- npointsAll[npointsAll[,2]==unique(npointsAll[,2])[k],]
      npointsk <- npoints[npoints[,2]==unique(npoints[,2])[k],]
      
      newAllk=matrix(nrow=sum(npointsAllk[,1]), ncol=3)
      newAllCnt=1
      newk=matrix(nrow=sum(npointsk[,1]), ncol=3)
      newCnt=1
      
      printTiming("Calc Probabilities", {
        for(i in 1:length(selectionGrid)){
          # probs for observation k of third variable from 3d density
          probsAll=MestimatesAll$estimate[cbind(match(grid[selectionGrid[[i]],1],
                                                      MestimatesAll$eval.points[[1]]),
                                                match(grid[selectionGrid[[i]],2],
                                                      MestimatesAll$eval.points[[2]]),
                                                match(MestimatesAll$eval.points[[3]][k],
                                                      MestimatesAll$eval.points[[3]]))]+1E-16
          probs=Mestimates$estimate[cbind(match(grid[selectionGrid[[i]],1],
                                                Mestimates$eval.points[[1]]),
                                          match(grid[selectionGrid[[i]],2],
                                                Mestimates$eval.points[[2]]),
                                          match(Mestimates$eval.points[[3]][k],
                                                Mestimates$eval.points[[3]]))]+1E-16
          probsAll = probsAll/sum(probsAll)
          probs = probs/sum(probs)
          probsPROP=(probs/probsAll)*(sum(datak[i,3])/sum(datak[i,4]))
          points=matrix(ncol=2,grid[selectionGrid[[i]],])
          points=cbind(points,as.integer(unique(data[,5])[k]))
          if(length(selectionGrid[[i]])==0){points <- matrix(ncol=2,shapefile@polygons[[i]]@labpt)
          probs=1
          probsAll=1}
          if(npointsAllk[,1][i]>0){
            sampleAll=sample(1:nrow(points),size=max(0,npointsAllk[,1][i],na.rm=T),replace=T,prob=probsAll)
            newAllk[newAllCnt:(newAllCnt+npointsAllk[,1][i]-1), ] = points[sampleAll,]
            newAllCnt = newAllCnt+npointsAllk[,1][i]
          }
          if(npointsk[,1][i]>0){
            if (length(sampleAll) == 1) {
              sampleProp = sampleAll
            } else {
              sampleProp <- sample(sampleAll,
                                   size = max(0, npointsk[,1][i], na.rm = T), prob = probsPROP[sampleAll])
            }
            newk[newCnt:(newCnt+npointsk[,1][i]-1), ] = points[sampleProp,]
            newCnt = newCnt+npointsk[,1][i]
          }
        }
      }) # end of i
      
      
      newAllk <- data.frame(newAllk) %>% group_by_all() %>% count %>% distinct %>% as.matrix() # slightly faster but requires dplyr
      # newAllk <- as.matrix(ddply(data.frame(newAllk),.(X1,X2,X3),nrow)) # uses plyr
      
      newAll <- rbind(newAll,newAllk)
      new <- rbind(new,newk)
      
    }
    
    #recompute H
    printTiming("Recompute H", {
      if(adaptive==FALSE){
        H <- ks::Hpi(x = new, binned=TRUE) * 2
      }
      if(adaptive==TRUE){
        H <- ks::Hpi(x = new, binned=TRUE)
        # adjusted for 3rd dimension
        H <- sqrt(sqrt(H[1,1]*H[2,2]*H[3,3]))
      }
    })
    
    wAll = newAll[,4]
    w = rep(1, nrow(new))
    
    if(boundary == TRUE){
      printTiming("Calc Weights", {
        if(fastWeights == FALSE || j <= ceiling(0.1*(burnin+samples))){
          weights <- calcWeights_fast( inside  = inside ,
                                       outside = outside,
                                       gridx   = gridx  ,
                                       gridy   = gridy  ,
                                       H       = H[1:2,1:2])
          
        }
      })
      
      printTiming("Match Weights", {
        new_2SingleID    = apply(new, 1, function(x) { x[1]*Max_2SingleID + x[2] } )
        newAll_2SingleID = apply(newAll, 1, function(x) { x[1]*Max_2SingleID + x[2] } )
        w <- weights[ fmatch(new_2SingleID, inside_2SingleID, nomatch = 1) ]
        wAll <- weights[ fmatch(newAll_2SingleID, inside_2SingleID, nomatch = 1) ] * wAll
      })
    }
    #recompute density
    printTiming("Recompute Density", {
      if(adaptive==FALSE){
        # added new weight variable and adjusted x, H, gridsize, xmin and xmax for third dimension
        MestimatesAll <- kde(newAll[,1:3], H = H,
                             gridsize=c(length(gridx),length(gridy),length(unique(data[,5]))),
                             bgridsize=c(length(gridx),length(gridy),length(unique(data[,5]))),
                             xmin=c(min(gridx),min(gridy),min(unique(data[,5]))),
                             xmax=c(max(gridx),max(gridy),max(unique(data[,5]))),
                             binned=TRUE, 
                             w = wAll/mean(wAll))
        
        Mestimates <- ks::kde(x=new, H=H,
                              gridsize=c(length(gridx),length(gridy),length(unique(data[,5]))),
                              bgridsize=c(length(gridx),length(gridy),length(unique(data[,5]))),
                              xmin=c(min(gridx),min(gridy),min(unique(data[,5]))),
                              xmax=c(max(gridx),max(gridy),max(unique(data[,5]))),
                              binned=TRUE, w = w/mean(w))
      }
      if(adaptive==TRUE){
        counts <- plyr::count(new)
        MestimatesAd <- sparr::bivariate.density(data=counts[,c(1:2)],pilotH=H,res=length(gridx),xrange=range(gridx),
                                                 yrange=range(gridy),adaptive=TRUE,
                                                 comment=FALSE, counts=counts[,4])
        Mestimates$estimate=MestimatesAd$Zm
      }
      #weird behaviour of ks::kde returning negative densities (probably due to numeric problems)
      Mestimates$estimate[Mestimates$estimate < 0] <- 0
      MestimatesAll$estimate[MestimatesAll$estimate < 0] <- 0
      
    })
    
    printTiming("Delete Shapes", {
      if(!is.null(deleteShapes)){
        for(k in 1:length(MestimatesAll$eval.points[[3]])){
          Mestimates$estimate[,,k][-setdiff(unlist(selectionGrid), unselectionGrid)] = 0
          MestimatesAll$estimate[,,k][-setdiff(unlist(selectionGrid), unselectionGrid)] = 0
        }
      } else{
        for(k in 1:length(MestimatesAll$eval.points[[3]])){
          Mestimates$estimate[,,k][-(unlist(selectionGrid))] = 0
          MestimatesAll$estimate[,,k][-(unlist(selectionGrid))] = 0
        }
      }
    })
    
    printTiming("Match Density", {
      # adjusted for 3rd dimension
      densityAll <- array(NA,dim(MestimatesAll$estimate))
      densityPart <- array(NA,dim(Mestimates$estimate))
      for(k in 1:length(MestimatesAll$eval.points[[3]])){
        densityAll[cbind(match(grid[selectionGrid %>% unlist ,1],MestimatesAll$eval.points[[1]]),
                         match(grid[selectionGrid %>% unlist,2],MestimatesAll$eval.points[[2]]),
                         match(MestimatesAll$eval.points[[3]][k],MestimatesAll$eval.points[[3]]))] =
          MestimatesAll$estimate[cbind(match(grid[selectionGrid %>% unlist ,1],MestimatesAll$eval.points[[1]]),
                                       match(grid[selectionGrid %>% unlist,2],MestimatesAll$eval.points[[2]]),
                                       match(MestimatesAll$eval.points[[3]][k],MestimatesAll$eval.points[[3]]))]
        densityPart[cbind(match(grid[selectionGrid %>% unlist ,1],MestimatesAll$eval.points[[1]]),
                          match(grid[selectionGrid %>% unlist,2],MestimatesAll$eval.points[[2]]),
                          match(MestimatesAll$eval.points[[3]][k],MestimatesAll$eval.points[[3]]))] =
          Mestimates$estimate[cbind(match(grid[selectionGrid %>% unlist ,1],MestimatesAll$eval.points[[1]]),
                                    match(grid[selectionGrid %>% unlist,2],MestimatesAll$eval.points[[2]]),
                                    match(MestimatesAll$eval.points[[3]][k],MestimatesAll$eval.points[[3]]))]
      }
    })
    
    printTiming("Assignments", {
      # adjusted for 3rd dimension
      proportion <- array(NA,dim(Mestimates$estimate))
      for(k in 1:length(MestimatesAll$eval.points[[3]])){
        densityPart[,,k] <- densityPart[,,k]/sum(densityPart[,,k], na.rm=TRUE)
        densityAll[,,k] <- densityAll[,,k]/sum(densityAll[,,k], na.rm=TRUE)
        densityAll[,,k] <- densityAll[,,k]+1E-96
        proportion[,,k] <- (densityPart[,,k]/densityAll[,,k]) * (sum(data[data[,5]==unique(data[,5])[k],][,3])/sum(data[data[,5]==unique(data[,5])[k],][,4]))
      }
      proportion[proportion<0]=0
      proportion[proportion>1]=1
      
      Mestimates$estimate[is.na(Mestimates$estimate)]=1E-96
      ret$resultDensity[j,,,] <- Mestimates$estimate
      ret$proportionArray[j,,,] <- proportion
      ret$resultX[j,,] <- new
    })
    printInfo("Iteration: ",j," of ", burnin+samples," in chain ", chain)
  }
  ret$Mestimates = Mestimates
  
  if (is.na(saveAsBaseFileName))
  {
    return( ret )
  }
  else
  {
    save(ret, file = paste(saveAsBaseFileName, chain, sep=""), envir = environment())
    return( NA )
  }
}


#' Transfer observations to other shape
#' @param Mestimates Estimation object created by functions dshapebivr and dbivr
#' @param shapefile The new shapefile for which the observations shall be transferred to
#' @return
#' The function returns the count, sd and 90%-coverage interval of the observations in the different shapfile:
#' @export
toOtherShape <- function(Mestimates, shapefile) {
  meanCount <- rep(NA, length(shapefile@polygons))
  sdCount <- rep(NA, length(shapefile@polygons))
  q0.05Count <- rep(NA, length(shapefile@polygons))
  q0.95Count <- rep(NA, length(shapefile@polygons))
  
  for (i in 1:length(shapefile@polygons)) {
    counts <-
      Reduce(c, lapply(1:dim(Mestimates$resultX)[1], function(j) {
        sapply((Mestimates$burnin + 1):(Mestimates$burnin + Mestimates$samples),
               function(x) {
                 lapply(1:length(shapefile@polygons[[i]]@Polygons), function(k) {
                   if (shapefile@polygons[[i]]@Polygons[[k]]@hole == FALSE) {
                     if(length(dim(Mestimates$resultX)) == 3){
                       which(
                         point.in.polygon(
                           Mestimates$resultX[x, , 1],
                           Mestimates$resultX[x, , 2],
                           shapefile@polygons[[i]]@Polygons[[k]]@coords[, 1],
                           shapefile@polygons[[i]]@Polygons[[k]]@coords[, 2]) == 1
                       )
                     } else {
                     which(
                       point.in.polygon(
                         Mestimates$resultX[j, x, , 1],
                         Mestimates$resultX[j, x, , 2],
                         shapefile@polygons[[i]]@Polygons[[k]]@coords[, 1],
                         shapefile@polygons[[i]]@Polygons[[k]]@coords[, 2]) == 1
                     )
                    }
                   }
                 }) %>% unlist %>% unique %>% length
               })
      }))
    
    meanCount[i] <- round(mean(counts))
    sdCount[i] <- round(sd(counts))
    q0.05Count[i] <- round(quantile(counts, 0.05))
    q0.95Count[i] <- round(quantile(counts, 0.95))
  }
  return(data.frame(meanCount, sdCount, q0.05Count, q0.95Count))
}

setKernelheapingLoglevel <- function(logInfo = TRUE, logDump = FALSE, logTiming = FALSE) {
  assign("logLevelInfo"  , logInfo  ,  envir=kernelheapingEnv)
  assign("logLevelDump"  , logDump  ,  envir=kernelheapingEnv)
  assign("logLevelTiming", logTiming,  envir=kernelheapingEnv)
}

setKernelheapingOutputLevel <- function(outputLevel) {
  assign("outputLevel"  , outputLevel  ,  envir=kernelheapingEnv)
}

isOutputlevelHigh <- function() {
  if ( exists("outputLevel", envir=kernelheapingEnv) && get("outputLevel", envir=kernelheapingEnv) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

printInfo <- function(...) { if ( exists("logLevelInfo", envir=kernelheapingEnv) && get("logLevelInfo", envir=kernelheapingEnv)) { print( paste(sep = "", ...)  ) } }
printDump <- function(...) { if ( exists("logLevelDump", envir=kernelheapingEnv) && get("logLevelDump", envir=kernelheapingEnv)) { print( paste(sep = "", ...)  ) } }

printTiming <- function(descr, ...) { 
  if ( exists("logLevelTiming", envir=kernelheapingEnv) && get("logLevelTiming", envir=kernelheapingEnv)) { 
    t = system.time( ... )
    print( paste(sep = "", "Timing ", descr, ": ", round(t[1], digits=2), " (", round(t[3], digits=2), ")")  ) 
  } else {
    (...)
  }
}

