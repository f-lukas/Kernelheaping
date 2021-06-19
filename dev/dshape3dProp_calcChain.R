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
  # adjusted for 3rd dimension
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
    
    # subset data and npoints, keeping only time period k
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
        # adjusted for 3rd dimension
        # probs for time period k from 3d density
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
    
    # binding sampled geocoordinates for each time period
    # BREAKING change: due to the third dimension matrix 'newAll' and 'new' grow with the factor equal to number of time periods
    # this can be problematic especially if number of observations for complete observations are large or the number of time periods are large
    # in case of the german population ~83 million and three time periods, matrix 'newAll' can grow to ~4gb or even more
    # therefore: keeping only unique rows (combinations of geocoordinates and time period) and adding a count column
    # later we use the count column as weight variable to estimate the new density
    # (currently) the count column is only added to matrix 'newAll'
    # matrix 'new' still contains duplicate rows (one row for every observation) as in the 2d case of dshapebivrProp()
    # while matrix 'new' is used to recompute bandwidth H, the method is equal/comparable to the 2 dimensional case of dshapebivrProp()
    
    # in case of boundary = T, weights are needed to adjust boundaries. At this time boundary weights will be multiplied with count weights.
    
    
    newAllk <- data.frame(newAllk) %>% group_by_all() %>% count %>% distinct %>% as.matrix() # slightly faster but requires dplyr
    # newAllk <- as.matrix(ddply(data.frame(newAllk),.(X1,X2,X3),nrow))
    
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