## Quickly get country centroids 
library(rgeos)
library(rworldmap)
# get world map
wmap <- getMap(resolution="high")

# get centroids
centroids <- gCentroid(wmap, byid=TRUE)

# get a data.frame with centroids
ccoords <- as.data.frame(centroids)
head(ccoords)

ccoords$country <- row.names(df)