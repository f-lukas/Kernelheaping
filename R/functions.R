#' Conditional Posterior for X given W, beta, a, p
#' @param rounds rounding values
#' @param Bias Bias parameter (on probit scale)
#' @param RR Threshold values for rounding parameters
#' @param beta acceleration parameter
#' @param gridx grid on which density is evaluated
#' @return List with Probabilities
#' @export
Rprx2=function(rounds, Bias, RR, beta, gridx){
  rprobs=rprobsRR2(RR,beta,gridx)
  b1=matrix(0,ncol=length(gridx),nrow=length(rounds))
  b2=matrix(0,ncol=length(gridx),nrow=length(rounds))
  possW=seq(from=plyr::round_any(min(gridx),min(rounds)),to=plyr::round_any(max(gridx),min(rounds)),by=min(rounds))
  numbers=matrix(0,nrow=length(rounds)*2,ncol=length(possW))
  for(i in 1:length(rounds)){
    b1[i,which(gridx%%rounds[i]<0.5*rounds[i])]=b1[i,which(gridx%%rounds[i]<0.5*rounds[i])]+Bias*rprobs[i,which(gridx%%rounds[i]<0.5*rounds[i])]/length(gridx)
    b2[i,which(gridx%%rounds[i]>0.5*rounds[i])]=b2[i,which(gridx%%rounds[i]>0.5*rounds[i])]+(1-Bias)*rprobs[i,which(gridx%%rounds[i]>0.5*rounds[i])]/length(gridx)
  }
  b=rbind(b1,b2)
  
  b=apply(b,2,function(x) x/sum(x))
  colnames(b)=gridx
  for(i in 1:length(c(rounds,rounds))){
    atemp=factor(plyr::round_any(gridx,accuracy=c(rounds,rounds)[i],f=round),levels=possW)
    temp=tapply(b[i,],INDEX=list(atemp),FUN=function(x) sum(x,na.rm=T))
    temp[is.na(temp)]=0
    numbers[i,]=temp
  }
  colnames(b)=gridx
  colnames(numbers)=possW
  splitted=apply(numbers,2,function(x) x/sum(x))
  return(list(splitted,b))
}
rprobsRR2=function(RR,beta,gridx){sapply(1:length(gridx), function(x) diff(c(0,pnorm(RR+beta*log(abs(gridx[x])+0.00001)),1)))}

logLik2=function(par, new, rguesstemp, unequal, setBias, rounds, gridx){
  RR=sort(par[1:(length(rounds)-1)])
  Bias=0
  if(setBias==TRUE) {Bias=par[length(rounds)]}
  if(unequal==FALSE){
    pgivenX=Rprx2(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=0, gridx=gridx)[[2]]
    return(-(sum(log(sapply(1:length(new),function(x) (pgivenX[rguesstemp[x],as.character(new[x])]))))))
  }
  if(unequal==TRUE) {
    beta=par[length(par)]
    pgivenX=Rprx2(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=beta, gridx=gridx)[[2]]
    return(-(sum(log(sapply(1:length(new),function(x) (pgivenX[rguesstemp[x],as.character(new[x])]))))))
  }
}

#' Kernel density estimation for heaped data
#' @param xheaped heaped values from which to estimate density of x
#' @param rounds rounding values 
#' @param burnin burn-in sample size
#' @param samples sampling iteration size
#' @param setBias if TRUE a rounding Bias parameter is estimated. For values above 0.5, the respondents
#' are more prone to round down, while for values < 0.5 they are more likely to round up  
#' @param bw bandwidth selector method, defaults to "nrd0" see \code{density} for more options  
#' @param boundary TRUE for positive only data (no positive density for negative values)
#' @param unequal if TRUE a probit model is fitted for the rounding probabilities with log(true value) as regressor
#' @param adjust as in \code{density}, the user can multiply the bandwidth by a certain factor such that bw=adjust*bw
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
#' ###################
#' #Simulate Data
#' ###################
#' Sim1 <- createSim.Kernelheaping(n=500, distribution="norm",rounds=c(1,10,100),
#' rprobs=c(0.3,0.4,0.3), sd=100)
#' \dontrun{est <- dheaping(Sim1$xheaped,rounds=Sim1$rounds)
#' plot(est,trueX=Sim1$x)}
#' 
#' #Biased rounding
#' Sim2 <- createSim.Kernelheaping(n=500, distribution="gamma",rounds=c(1,2,5,10),
#'                      rprobs=c(0.1,0.15,0.4,0.35),downbias=0.2,shape=4,scale=8,offset=45)
#' \dontrun{est <- dheaping(Sim2$xheaped, rounds=Sim2$rounds, setBias=T, bw="SJ")
#' plot(est, trueX=Sim2$x)
#' summary(est)
#' tracePlots(est)}
#' 
#' Sim3 <- createSim.Kernelheaping(n=500, distribution="gamma",rounds=c(1,2,5,10),
#' rprobs=c(0.8,0.1,0.05,0.05), downbias=0.75, Beta=-1, shape=4, scale=8)
#' \dontrun{est <- dheaping(Sim3$xheaped,rounds=Sim3$rounds,boundary=TRUE,unequal=TRUE,setBias=T)
#' plot(est,trueX=Sim3$x)}
#' 
#' ###################
#' #Real Data Example
#' ###################
#' # Student learning hours per week
#' data(students)
#' xheaped <- as.numeric(na.omit(students$StudyHrs))
#' \dontrun{est <- dheaping(xheaped,rounds=c(1,2,5,10), boundary=TRUE, unequal=TRUE)
#' plot(est)
#' summary(est)}
#' @export
dheaping <- function(xheaped, rounds, burnin=5, samples=10, setBias=FALSE,
                             bw= "nrd0", boundary=FALSE, unequal=FALSE, adjust=1){
  xheaped=plyr::round_any(xheaped,accuracy=min(rounds))
  #Create grid
  gridx=seq(min(xheaped)-max(rounds)*0.5+0.05*min(rounds),max(xheaped)+max(rounds)*0.55-0.05*min(rounds),by=min(rounds)*0.1)
  #Result matrices
  resultDensity=matrix(nrow=samples+burnin,ncol=length(gridx))
  resultX=matrix(nrow=samples+burnin,ncol=length(xheaped))
  resultRR=matrix(nrow=samples+burnin,ncol=length(rounds)-1)
  resultBias=matrix(nrow=samples+burnin,ncol=1)
  resultBeta=matrix(nrow=samples+burnin,ncol=1)
  #Pilot Estimation
  if(boundary==FALSE){Mestimates <- density(xheaped,from=min(gridx),to=max(gridx),n=length(gridx),bw=max(rounds)/4)$y}
  if(boundary==TRUE){Mestimates <- evmix::dbckden(gridx,xheaped,bw=max(rounds)/4,bcmethod="simple")}
  #Starting values for x(new), Rounding values(rguesstemp), und Rounding thresholds(RR), Bias (Bias) and beta
  rguess=sapply(1:length(xheaped),function(x) rounds[max(which(round(xheaped[x])%%round(rounds)==0))])
  rguessx=factor(xheaped%%max(rounds),levels=seq(from=0,to=max(rounds)-min(rounds),by=min(rounds)))
  rguesstemp=rguess
  new=xheaped
  RR=qnorm(seq(1:(length(rounds)-1))/(length(rounds)))
  Bias=0
  beta=0
  for(j in 1:(burnin+samples)){
    Rprs=Rprx2(rounds=rounds,RR=RR, Bias=pnorm(Bias), beta=beta, gridx=gridx)
    for(i in 1:length(xheaped)){
      selectionupper=lapply(1:length(rounds), function(x) 
        which(gridx>xheaped[i]-rounds[x]*0.5&gridx<=xheaped[i]))
      selectionlower=lapply(1:length(rounds), function(x)
        which(gridx<xheaped[i]+rounds[x]*0.5&gridx>=xheaped[i]))
      selection=c(selectionlower,selectionupper)
      selectionprobs=lapply(1:length(selection),function(x) 
        Mestimates[selection[[x]]]*
          Rprs[[2]][x,as.character(gridx[selection[[x]]])]/
          (0.00000001+sum(Rprs[[2]][x,as.character(gridx[selection[[x]]])]))*
          (Rprs[[1]][x,as.character(xheaped[i])]))
      listlength=lapply(selection, function(x) length(x))
      temprounds=unlist(lapply(1:length(selection), function(x) rep(x,times=listlength[[x]])))
      temp=sample(1:length(unlist(selectionprobs)),size=1,prob=unlist(selectionprobs))
      rguesstemp[i]=temprounds[temp]
      new[i]=gridx[unlist(selection)[temp]]
    }
    #Laplace-Approximation:
    par=RR
    if(setBias==TRUE){par=c(par,Bias)}
    if(unequal==TRUE){par=c(par,beta)}
    laplace=optim(par,logLik2,new=new,rguesstemp=rguesstemp, unequal=unequal,
                  setBias=setBias, rounds=rounds, gridx=gridx ,hessian=T,method="BFGS")
    par=MASS::mvrnorm(1,mu=laplace$par,Sigma=solve(laplace$hessian+diag(0.0001,length(par))))
    if(setBias==TRUE) {Bias=par[length(rounds)]}
    if(unequal==TRUE) {beta=par[length(par)]}
    RR=sort(par[1:(length(rounds)-1)])
    h <- density(new,bw=bw)$bw*adjust
    if(boundary==TRUE){Mestimates <- evmix::dbckden(gridx,new,bw=h,bcmethod="simple")}
    if(boundary==FALSE){Mestimates <- density(new,from=min(gridx),to=max(gridx),n=length(gridx),bw=h)$y}
    #Save Results
    resultDensity[j,]=Mestimates
    resultRR[j,]=RR
    resultBias[j,]=pnorm(Bias)
    resultBeta[j,]=beta
    resultX[j,]=new
    print(paste("Iteration:",j,"of", burnin+samples))
  }
  meanPostDensity=apply(resultDensity[-c(1:burnin),],2,mean)
  est<-list(meanPostDensity=meanPostDensity,resultDensity=resultDensity,resultRR=resultRR, resultBias=resultBias,
            resultBeta=resultBeta, resultX=resultX, xheaped=xheaped, gridx=gridx, boundary=boundary, rounds=rounds, setBias=setBias, unequal=unequal, bw=bw, burnin=burnin, samples=samples)
  class(est) <- "Kernelheaping"
  return(est)
}

#' Plot Kernel density estimate of heaped data naively and corrected by partly bayesian model
#' @param x Kernelheaping object produced by \code{dheaping} function
#' @param trueX optional, if true values X are known (in simulations, for example) the 'Oracle' density estimate is added as well
#' @param ... additional arguments given to standard plot function
#' @return plot with Kernel density estimates (Naive, Corrected and True (if provided))
#' @export
plot.Kernelheaping <- function(x,trueX=NULL, ...){
  if(x$boundary==FALSE){plot(density(x$xheaped,bw=x$bw),xlab="x",ylab="Density", ...)}
  if(x$boundary==TRUE){plot(evmix::dbckden(x$gridx,x$xheaped,bw=density(x$xheaped,bw=x$bw)$bw,bcmethod="simple")~x$gridx,type="l",xlab="x",ylab="Density", ...)}
  lines(x$meanPostDensity~x$gridx,col="red")
  if(!is.null(trueX)){lines(density(trueX),col="blue")}
  if(is.null(trueX)){legend("topright",c("Naive","Corrected"),col=c("black","red"),lty=1)}
  if(!is.null(trueX)){legend("topright",c("Naive","Corrected","Oracle"),col=c("black","red","blue"),lty=1)}
}

#' Prints some descriptive statistics (means and quantiles) for the estimated rounding, bias and acceleration (beta) parameters
#' @param object Kernelheaping object produced by \code{dheaping} function
#' @param ... unused
#' @return Prints summary statistics
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
  par(mfrow=c(1,ncol(x$resultRR)))
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
}

#' Simulation of heaping correction method
#' @param simRuns number of simulations runs
#' @param n sample size
#' @param distribution name of the distribution where random sampling is available, e.g. "norm"
#' @param rounds rounding values
#' @param rprobs rounding probabilities (for Beta=0)
#' @param ... additional attributes handed over to \code{createSim.Kernelheaping} or \code{dheaping}
#' @return List of estimation results
#' @examples
#' \dontrun{Sims1 <- sim.Kernelheaping(simRuns=2, n=500, distribution="norm", 
#' rounds=c(1,10,100), rprobs=c(0.3,0.4,0.3), sd=100)}
#' @export
sim.Kernelheaping <- function(simRuns, n, distribution, rounds, rprobs, ...){
  lapply(1:simRuns, function(x){
    print(x)
    Sim=createSim.Kernelheaping(n, distribution, rounds, rprobs, ...)
    est=dheaping(Sim$xheaped, rounds=Sim$rounds)
    return(list(est,Sim))
  })
}

#' Create heaped data for Simulation
#' @param n sample size
#' @param distribution name of the distribution where random sampling is available, e.g. "norm"
#' @param rounds rounding values
#' @param rprobs rounding probabilities (for Beta=0)
#' @param offset certain value added to all observed random samples
#' @param downbias bias parameter
#' @param Beta acceleration paramter
#' @param ... additional attributes handed over to "rdistribution" (i.e. rnorm, rgamma,..)
#' @return List of heaped values, true values and input parameters
#' @export
createSim.Kernelheaping <- function(n, distribution, rounds, rprobs, offset=0, downbias=0.5, Beta=0, ...){
  x=do.call(paste("r",distribution,sep=""), list(n, ...))
  RR0=qnorm(cumsum(rprobs))[-length(rprobs)]
  down=lapply(x,function(x) (x%%rounds<=rounds*0.5)*downbias*rprobsRR2(RR0,Beta,x)+(x%%rounds>=rounds*0.5)*(1-downbias)*rprobsRR2(RR0,Beta,x))
  down=lapply(down, function(x) x/sum(x))
  roundings=sapply(down,function(z) sample(rounds,size=1,replace=TRUE,prob=z))
  xheaped=plyr::round_any(x,accuracy=roundings, f=round)
  return(list(xheaped=xheaped, rounds=rounds, rprobs=rprobs, downbias=downbias, Beta=Beta, x=x))
}

