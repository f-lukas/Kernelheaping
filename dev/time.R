
library(Kernelheaping)
# library(MASS)
# library(ks)
# library(sparr)
library(sp)
library(plyr)
library(fastmatch)
library(magrittr)
library(mvtnorm)
library(parallel)
library(dplyr)


# set wd to path of current file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load data
# data.RData contains
# shapefile
# data_3days: data.frame with data for 15.03.2021 - 17.03.2021
# data_7days: data.frame with data for 15.03.2021 - 21.03.2021
# data_14days: data.frame with data for 15.03.2021 - 28.03.2021
# data_28days: data.frame with data for 01.03.2021 - 28.03.2021
load(paste0(getwd(),"/data.RData"))

# load kernelheaping functions
source(paste0(dirname(getwd()),"/R/functions.R"))
source(paste0(getwd(),"/dshape3dProp.R"))
source(paste0(getwd(),"/dshape3dProp_calcChain.R"))

######################################################

# create dataframe for duration
timeDf<-NULL

######################################################

# create function
saveTime <- function(data,burnin,samples,numChains,numThreads){

start <- Sys.time()
set.seed(1)
est <- dshape3dProp(data = data,
                    burnin = burnin,
                    samples = samples,
                    shapefile = shapefile,
                    gridsize = 200,
                    numChains = numChains,
                    numThreads = numThreads)
end <- Sys.time()
print(end-start)
timeDf <- timeDf %>% bind_rows(data.frame(nDates = length(unique(data$date)),
                                          time = as.numeric(end-start),
                                          burnin = burnin,
                                          samples = samples,
                                          chains = numChains,
                                          threads = numThreads))
return(timeDf)
}

######################################################

timeDf <- saveTime(data = data_7days, burnin = 3, samples = 7, numChains = 1, numThreads = 1)
timeDf <- saveTime(data = data_14days, burnin = 3, samples = 7, numChains = 1, numThreads = 1)
timeDf <- saveTime(data = data_28days, burnin = 3, samples = 7, numChains = 1, numThreads = 1)
timeDf <- saveTime(data = data_28days, burnin = 5, samples = 15, numChains = 1, numThreads = 1)

timeDf <- saveTime(data = data_14days, burnin = 3, samples = 7, numChains = 16, numThreads = 1)
timeDf <- saveTime(data = data_14days, burnin = 3, samples = 7, numChains = 16, numThreads = 8)
timeDf <- saveTime(data = data_14days, burnin = 5, samples = 15, numChains = 16, numThreads = 8)


write.csv(timeDf,"timeDf.csv",row.names = F)


