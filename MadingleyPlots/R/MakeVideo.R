
MakeVideo <- function(resultsDir,videoName,outDir,
                      label=NULL,
                      vars=c("Biomass density"),
                      map="World",frameRate=10){
  
  if (suppressWarnings(system(
    "ffmpeg",show.output.on.console = FALSE)==127)){
    stop("You need to install ffmpeg")
  }
  
  if (map == "Africa"){
  } else {
    data(World)
  }
  
  
  LogVariables<-list(
    "Abundance density" = TRUE,
    "Biomass density" = TRUE,
    "autotrophbiomass density" = TRUE,
    "carnivoreabundance density" = TRUE,
    "carnivorebiomass density" = TRUE,
    "herbivoreabundance density" = TRUE,
    "herbivorebiomass density" = TRUE,
    "omnivoreabundance density" = TRUE,
    "omnivorebiomass density" = TRUE
    
  )
  
  GramsBiomassVariables<-list(
    "Abundance density" = FALSE,
    "Biomass density" = TRUE,
    "autotrophbiomass density" = TRUE,
    "carnivoreabundance density" = FALSE,
    "carnivorebiomass density" = TRUE,
    "herbivoreabundance density" = FALSE,
    "herbivorebiomass density" = TRUE,
    "omnivoreabundance density" = FALSE,
    "omnivorebiomass density" = TRUE
    
  )
  
  PlantBiomassVariables<-list(
    "Abundance density" = FALSE,
    "Biomass density" = FALSE,
    "autotrophbiomass density" = TRUE,
    "carnivoreabundance density" = FALSE,
    "carnivorebiomass density" = FALSE,
    "herbivoreabundance density" = FALSE,
    "herbivorebiomass density" = FALSE,
    "omnivoreabundance density" = FALSE,
    "omnivorebiomass density" = FALSE
    
  )
  
  stopifnot(all(vars %in% names(LogVariables)))
  stopifnot(all(vars %in% names(GramsBiomassVariables)))
  stopifnot(all(vars %in% names(PlantBiomassVariables)))
  
  varsToTransform<-which(unlist(LogVariables[match(vars,names(LogVariables))]))
  varsToTransformKG<-which(unlist(GramsBiomassVariables[
    match(vars,names(GramsBiomassVariables))]))
  varsToTransformPlantBiomass<-which(unlist(PlantBiomassVariables[match(
    vars,names(PlantBiomassVariables))]))
  
  .Log("Finding Madingley output files\n")
  files<-.ListGridOutputFiles(resultsDir)
  
  # Find the simulation numbers
  sims.re<-regexpr("_[0-9]+.nc",files)
  sims<-as.list(unique(substr(files,sims.re,sims.re+
                                attr(sims.re,"match.length")-1)))
  
  if(is.null(label)){
    label<-unique(substr(files,1,sims.re-1))
    stopifnot(length(label)==1)
    label<-label[1]
  } else {
    label <- paste("GridOutputs_",label,sep="")
  }
  
  .Log(paste("Found results for ",length(sims)," simulations\n",sep=""))
  
  .Log("Getting basic information about simulations\n")
  sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],sep="")
  data<-open.sds(sds.path)
  long<-get.sds(data,"Longitude")
  lat<-get.sds(data,"Latitude")
  ar<-length(lat)/length(long)
  
  allTimes <- get.sds(data,"Time step")
  
  if (dir.exists(paste(outDir,"/frames_",videoName,sep=""))){
    unlink(paste(outDir,"/frames_",videoName,sep=""))
  }
  
  dir.create(paste(outDir,"/frames_",videoName,sep=""))
  
  .Log("Processing results\n")
  
  # Create matrices to hold the results for each specified variable
  allResults<-array(data = NA,dim = c(length(sims),length(allTimes),
                                      length(vars),length(long),length(lat)))
  
  # Loop over simulations in the ensemble
  s<-1
  for (sim in sims){
    sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,sep="")
    
    data<-open.sds(sds.path)
    
    # Populate the results matrices
    v<-1
    for (var in vars){
      allResults[s,,v,,]<-get.sds(data,var)
      v<-v+1
    }
    s<-s+1
  }
  
  allResults[is.na(allResults)]<-0
  
  # Calcualte the mean value across simulations
  resultsSimMean<-apply(allResults,c(2,3,4,5),mean,na.rm=TRUE)
  
  # Transform any variables stored in log space
  resultsSimMean[,varsToTransform,,]<-exp(resultsSimMean[,varsToTransform,,])-1
  
  # Convert any variables stored in grams to kilograms
  resultsSimMean[,varsToTransformKG,,]<-resultsSimMean[,varsToTransformKG,,]/1000
  
  # Convert any plant biomass variables to total biomass
  # (not just leaf biomass)
  for (i in varsToTransformPlantBiomass){
    for (s in 1:length(sims)){
      resultsSimMean[s,i,,]<-resultsSimMean[s,i,,]/
        .GetLPRatios(resultsDir)
    }
    
  }
  
  # Sum values across the specified variables
  resultsTotal<-apply(resultsSimMean,c(1,3,4),sum,na.rm=TRUE)
  
  resultsTotal[resultsTotal==0]<-NA
  
  # Base colour breaks on last time step
  brks <- quantile(resultsTotal[length(allTimes),,],probs = seq(from=0,to=1,length.out=11),
                   na.rm = TRUE)
  brks[1] <- -9e99
  brks[11] <- 9e99
  
  .Log("Plotting frames\n")
  for (t in 1:length(allTimes)){
    png(filename = paste(outDir,"/frames_",videoName,"/time",sprintf("%04d",t),
                         ".png",sep=""),
        width = dim(resultsTotal)[2]*10,height=dim(resultsTotal)[3]*10,
        units="px")
    
    resultsMap<-SpatialGridDataFrame(grid = GridTopology(
      cellcentre.offset = c(min(long),min(lat)),
      cellsize = c(unique(diff(long)),unique(diff(lat))),
      cells.dim = c(length(long),length(lat))),
      data = data.frame(band1=as.vector(resultsTotal[t,,length(lat):1])))
    
    image(resultsMap,breaks = brks,col=rev(brewer.pal(n = 10,name = "RdYlBu")),
          xaxt="n",yaxt="n",bty="n")
    
    if (map=="Africa"){
      AfricaOutline(add = TRUE)
    } else {
      plot(outline,add=TRUE)
    }
    
    
    dev.off()
  }
  
  .Log("Making Video\n")
  system(paste("ffmpeg -framerate ",frameRate," -i ",getwd(),"/",outDir,
               "frames_",videoName,"/time%04d.png -r 30 -vcodec mpeg1video ",getwd(),"/",outDir,videoName,
               ".mp4",sep=""))
  
}