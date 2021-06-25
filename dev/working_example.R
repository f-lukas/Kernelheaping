
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


data <- data_3days
head(data)
# convert dates back to date class (dshape3dProp needs dates as numeric)
# as.Date(unique(data$date), origin = "1970-01-01")

# not yet implemented or tested:
  # adaptive = TRUE
  # delete shapes

  set.seed(1)
  est <- dshape3dProp(data = data,
                          burnin = 2,
                          samples = 5,
                          shapefile = shapefile,
                          gridsize = 200,
                          numChains = 1,
                          numThreads = 1)


######################################################
# create map

library(ggplot2)

dateIndex <- 1 # adjust for different days

map.data <- data.frame(expand.grid(long = est$Mestimates$eval.points[[1]],
                                   lat = est$Mestimates$eval.points[[2]]),
                                   Density = est$proportions[,,dateIndex] %>% as.vector * 100000,
                                   date = as.Date(est$Mestimates$eval.points[[3]][dateIndex], origin = "1970-01-01")) %>% filter(Density > 0)


breaks <- c(0,35,50,100,150,200,300,500,2000)
nbreaks <- length(breaks)-1
cuts = format(breaks,big.mark=".",decimal.mark=",",digits=0,scientific = FALSE)
labels= c(paste0(" ≤",cuts[-1][1:nbreaks-1]),paste0(" >",cuts[-1][nbreaks-1]))
# based on RdYlGn palette (RColorBrewer) but with one color changed
cols <- c("#1A9850","#66BD63","#A6D96A","#EFF8B1","#FEE08B","#FDAE61","#F46D43","#D73027")

map.data$wert_breaks <- cut(map.data$Density,
                                breaks = breaks,
                                labels = labels,
                                right=FALSE,
                                include.lowest = T)


p <- ggplot(map.data) +
  coord_fixed(ratio=(1/cos(pi*51.163375/180)))+
  geom_raster(aes(long, lat, fill = wert_breaks)) +
  ggtitle(paste0("Lokale 7-Tage-Inzidenz\n",format(unique(map.data$date)-6,"%d.%m.%Y")," bis ",format(unique(map.data$date),"%d.%m.%Y"))) +
  labs(caption=expression('© GeoBasis-DE / BKG 2021 (Daten verändert); Robert Koch-Institut (RKI), dl-de/by-2-0;\nempirica regio (© Statistische Ämter des Bundes und der Länder, Deutschland, 2018-2021, dl-de/by-2-0, https://www.govdata.de/dl-de/by-2-0)'))+
  scale_fill_manual(
    values = cols,
    name = "Inzidenz (Fälle je 100.000 Einwohner)",
    drop = FALSE,
    guide = guide_legend(
      direction = "horizontal",
      keyheight = unit(4, units = "mm"),
      keywidth = unit(166 / length(labels), units = "mm"),
      title.position = 'top',
      title.hjust = 0.5,
      label.hjust = 0.5,
      nrow = 1,
      byrow = T,
      label.position = "bottom"))+
  theme(panel.grid.major = element_line(colour = "transparent"),
        rect = element_blank(),
        plot.title = element_text(size = 30, face = "bold"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position="bottom",
        plot.caption=element_text(vjust = -5),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        plot.margin=margin(t=20,b=20))

p
        