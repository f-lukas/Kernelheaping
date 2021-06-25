dshape3dProp <- function(data, burnin=2, samples=5, adaptive=FALSE,
                           shapefile, gridsize=200, boundary = FALSE,
                           deleteShapes = NULL, fastWeights = TRUE,
                           numChains=1, numThreads=1){
  ###########################################################
  ##### get polygon shape coordinates from data #####
  ###########################################################
  
  # order data
  data <- data[order(data[,5]),]
  
  pol.x <- list()
  pol.y <- list()
  for (i in 1:length(shapefile@polygons)) {
    pol.x[[i]] <- shapefile@polygons[[i]]@Polygons[[1]]@coords[,1]
    pol.y[[i]] <- shapefile@polygons[[i]]@Polygons[[1]]@coords[,2]
  }
  
  npointsAll <- data[ ,c(4,5)] # number of observations for complete observations and time periods
  npoints <- data[ ,c(3,5)] # number of observations in area for partial population and time periods
  
  if (length(gridsize) == 1) { gridsize = c(gridsize, gridsize) }
  #Create grid - sparr package requires identical grid length on both axis:
  gridx=seq(shapefile@bbox[1,1],shapefile@bbox[1,2],length=gridsize[1])
  gridy=seq(shapefile@bbox[2,1],shapefile@bbox[2,2],length=gridsize[2])
  grid <- as.matrix(expand.grid(gridx,gridy))
  
  # calculate weights for for complete observations and each time period
  dataAllWeight <- sapply(unique(data[,5]), function(x){
    d <- data[data[,5]==x,]
    nrow(d)*d[,4]/sum(d[,4])
  })
  
  # calculate weights for for partial population and each time period
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
  # adjusted for 3rd dimension
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
      # adjusted for 3rd dimension
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
    clusterEvalQ(cl, { library(ks); library(mvtnorm); library(dplyr); library(fastmatch); library(Kernelheaping); source(paste0(dirname(getwd()),"/R/functions.R")) })
    
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
      
      # adjusted for 3rd dimension
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
  # adjusted for 3rd dimension
  indexGridColumns <- c( length( dim(resultDensity[,-c(1:burnin),,,]) ) -2,
                         length( dim(resultDensity[,-c(1:burnin),,,]) ) -1,
                         length( dim(resultDensity[,-c(1:burnin),,,]) ))
  Mestimates$estimate=apply(resultDensity[,-c(1:burnin),,,],indexGridColumns,mean)
  est<-list(Mestimates=Mestimates,resultDensity=resultDensity,
            data=data, gridx=gridx, gridy=gridy, shapefile=shapefile,
            burnin=burnin, samples=samples, adaptive=adaptive, numChains=numChains,
            proportions = apply(proportionArray[,-c(1:burnin),,,],indexGridColumns,mean))
  if (isOutputlevelHigh()) { est[["resultX"]]=resultX }
  class(est) <- "bivshape"
  return(est)
}
