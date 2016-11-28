.Log <- function(...) {
  if(getOption('MPverbose', TRUE)) cat(...)
}

.ListCellOuputFiles<-function(resultsDir){
  
  files<-dir(resultsDir)
  files<-files[grep("BasicOutputs",files)]
  files<-files[grep("Cell",files)]
  
  return(files)
  
}

.ListGridOutputFiles<-function(resultsDir){
  
  files<-dir(resultsDir)
  files<-files[grep("GridOutputs",files)]
  
  return(files)
}

.ListMassBinsFiles<-function(resultsDir){
  
  files<-dir(resultsDir)
  files<-files[grep("MassBinsOutputs",files)]
  files<-files[grep("Cell",files)]
  
  return(files)
  
}

DegreeCellAreaKM <- function(lat, height, width) {
  # Returns the area in km squared of a grid cell in degrees of arc
  # lat - the latitudinal centre of the cell
  # height, width - the size of the grid cell in degrees
  
  # TODO Unit test
  # TODO Reference for this method
  
  radians <- function(theta) theta*pi/180.0
  
  # Convert the latitude into radians
  lat.rad <- radians(lat)
  
  # The equatorial and polar radii of the Earth in km
  eq.radius <-  6378137
  pol.radius <- 6356752.3142
  
  # Calculate cell area
  angular.eccentricity <- acos(radians(pol.radius/eq.radius))
  ecc.sq <- sin(radians(angular.eccentricity))^2
  flattening <- 1-cos(radians(angular.eccentricity))
  temp.val <- (eq.radius*cos(lat.rad))^2+(pol.radius*sin(lat.rad))^2
  m.phi <- ((eq.radius*pol.radius)^2)/(temp.val^1.5)
  n.phi <- (eq.radius^2)/sqrt(temp.val)
  lat.length <- pi/180*m.phi/1000
  long.length <- pi/180*cos(lat.rad)*n.phi/1000
  return (lat.length*height*long.length*width)
}


.GetLPRatios <- function(resultsDir){
  
  data(LeafToPlantRatio)
  
  if ("SpecificLocations.csv" %in% dir(resultsDir)){
    locations <- read.csv(paste(resultsDir,"/SpecificLocations.csv",sep=""))
    values <- extract(raster(lp.ratio),locations[,c(2,1)])
    
  } else {
    oldStyleInputs <- FALSE
    if ("SimulationControlParameters.csv" %in% dir(resultsDir)){
      initialization <- read.csv(paste(resultsDir,"/SimulationControlParameters.csv",sep=""))
    } else {
      initialization <- read.csv(paste(resultsDir,"/EcosystemModelInitialisation.csv",sep=""))
      oldStyleInputs <- TRUE
    }
    if (oldStyleInputs){
      ll <- as.numeric(paste(initialization$Value[(initialization$Parameter=="Leftmost Latitude")]))
      rl <- as.numeric(paste(initialization$Value[(initialization$Parameter=="Rightmost Latitude")]))
    } else {
      ll <- as.numeric(paste(initialization$Value[(initialization$Parameter=="Leftmost Longitude")]))
      rl <- as.numeric(paste(initialization$Value[(initialization$Parameter=="Rightmost Longitude")]))
    }
    bl <- as.numeric(paste(initialization$Value[(initialization$Parameter=="Bottom Latitude")]))
    tl <- as.numeric(paste(initialization$Value[(initialization$Parameter=="Top Latitude")]))
    cs <- as.numeric(paste(initialization$Value[(initialization$Parameter=="Grid Cell Size")]))
    
    print(ll)
    print(rl)
    print(bl)
    print(tl)
    print(cs)
    
    tempRaster <- raster(xmn=ll,xmx=rl,ymn=bl,ymx=tl,nrows=(tl-bl)/cs,ncols=(rl-ll)/cs)
    tempRaster <- resample(raster(lp.ratio),tempRaster,method="bilinear")
    values <- matrix(tempRaster@data@values,nrow = (rl-ll)/cs,ncol = (tl-bl)/cs)
    values <- values[,ncol(values):1]
  }
  
  values[is.na(values)]<-1
  
  return(values)
  rm(lp.ratio)
  
}
