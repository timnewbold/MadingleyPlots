PlotCells<-function(resultsDir,outDir=NULL,map="World"){
  
  if (map == "Africa"){
    data(Africa)
  } else {
    data(World)
  }
  
  locations<-read.csv(paste(resultsDir,"/SpecificLocations.csv",sep=""))
  
  if (!is.null(outDir)){
    pdf(paste(outDir,"CellsMap.pdf",sep=""),
        width = 12.5/2.54,height=6.25/2.54)
  }
    
  
  par(mar=c(0,0,0,0))
  
  plot(outline)
  
  points(locations$Longitude,locations$Latitude,pch=16,col="red")
  
  text(locations$Longitude,locations$Latitude,1:dim(locations)[1],pos=2)
  
  if (!is.null(outDir)) invisible(dev.off())
  
}