
path <- paste("msds:nc?file=","LeafToPlantBiomassRatios.nc",
              "&openMode=open",sep="")

data <- open.sds(path)

lp.ratio.grid <- get.sds(data,"SoilCMed")
lp.ratio.lons<-get.sds(data,"longitude")
lp.ratio.lats<-get.sds(data,"latitude")

lp.ratio <- SpatialGridDataFrame(grid = GridTopology(cellcentre.offset = c(
  min(lp.ratio.lons),min(lp.ratio.lats)),cellsize = c(unique(
    diff(lp.ratio.lons))[1],unique(diff(lp.ratio.lats))[1]),
  cells.dim = c(length(lp.ratio.lons),length(lp.ratio.lats))),
  data = data.frame(band1=as.vector(t(lp.ratio.grid[nrow(lp.ratio.grid):1,]))))

rm(path,data,lp.ratio.grid,lp.ratio.lats,lp.ratio.lons)
