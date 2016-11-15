
PlotGridCompare <- function(resultsDirBaseline,resultsDir,plotName,outDir,
                            labelBaseline=NULL,label=NULL,
                     vars=c("Biomass density"),
                     endTimeStepBaseline=NULL,numTimeStepsBaseline=12,
                     endTimeStep=NULL,numTimeSteps=12,
                     width=12.5,res=300,map="World",
                     returnMap=FALSE,
                     brks=NULL,col=NULL){
  
  if (map == "Africa"){
    data(Africa)
  } else {
    data(World)
  }
  
  if (any(!is.null(brks),!is.null(col))){
    if(!all(!is.null(brks),!is.null(col))){
      stop("Error: if one of brks or col is set, both must be")
    }
    if (length(brks) != length(col)+1){
      stop("Error: brks must have one more value than col")
    }
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
  filesBaseline<-.ListGridOutputFiles(resultsDirBaseline)
  
  # Find the simulation numbers
  sims.reBaseline<-regexpr("_[0-9]+.nc",filesBaseline)
  simsBaseline<-as.list(unique(substr(filesBaseline,sims.reBaseline,sims.reBaseline+
                                attr(sims.reBaseline,"match.length")-1)))
  
  if(is.null(labelBaseline)){
    labelBaseline<-unique(substr(filesBaseline,1,sims.reBaseline-1))
    stopifnot(length(labelBaseline)==1)
    labelBaseline<-labelBaseline[1]
  } else {
    labelBaseline <- paste("GridOutputs_",labelBaseline,sep="")
  }
  
  .Log(paste("Found baseline results for ",length(simsBaseline)," simulations\n",sep=""))
  
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
  
  .Log(paste("Found impact results for ",length(sims)," simulations\n",sep=""))
  
  if (length(sims) != length(simsBaseline)){
    .Log("Warning: number of simulations under impact not equal to number in baseline\n")
  }
  
  .Log("Getting basic information about simulations\n")
  sds.path<-paste("msds:nc?file=",resultsDirBaseline,"/",labelBaseline,simsBaseline[1],sep="")
  data<-open.sds(sds.path)
  longBaseline<-get.sds(data,"Longitude")
  latBaseline<-get.sds(data,"Latitude")
  
  allTimesBaseline <- get.sds(data,"Time step")
  if (is.null(endTimeStepBaseline)) endTimeStepBaseline <- tail(allTimesBaseline,1)
  timesBaseline <- (endTimeStepBaseline-numTimeStepsBaseline+1):endTimeStepBaseline
  
  sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],sep="")
  data<-open.sds(sds.path)
  long<-get.sds(data,"Longitude")
  lat<-get.sds(data,"Latitude")
  
  allTimes <- get.sds(data,"Time step")
  if (is.null(endTimeStep)) endTimeStep <- tail(allTimes,1)
  times <- (endTimeStep-numTimeSteps+1):endTimeStep
  
  stopifnot(all.equal(longBaseline,long))
  stopifnot(all.equal(latBaseline,lat))
  
  ar<-length(lat)/length(long)
  
  .Log("Initializing plot\n")
  tiff(paste(outDir,plotName,".tif",sep=""),width = width,height = width*ar,
       units = "cm",res=res,compression = "lzw")
  
  par(mar=c(0,0,0,0))
  
  .Log("Processing baseline results\n")
  
  # Create matrices to hold the results for each specified variable
  allResultsBaseline<-array(data = NA,dim = c(length(simsBaseline),length(allTimesBaseline),
                                      length(vars),length(longBaseline),length(latBaseline)))
  
  # Loop over simulations in the ensemble
  s<-1
  for (sim in simsBaseline){
    sds.path<-paste("msds:nc?file=",resultsDirBaseline,"/",labelBaseline,sim,sep="")
    
    data<-open.sds(sds.path)
    
    # Populate the results matrices
    v<-1
    for (var in vars){
      allResultsBaseline[s,,v,,]<-get.sds(data,var)
      v<-v+1
    }
    s<-s+1
  }
  
  # Select only the designated time steps
  allResultsBaseline<-allResultsBaseline[,match(timesBaseline,allTimesBaseline),,,,drop=FALSE]
  
  allResultsBaseline[is.na(allResultsBaseline)]<-0
  
  # Calcualte the mean value across these time steps for 
  # each variable within each simulation
  resultsTimesMeanBaseline<-apply(allResultsBaseline,c(1,3,4,5),mean,na.rm=TRUE)
#   # Calcualte the mean value across simulations
#   resultsSimMeanBaseline<-apply(resultsTimesMeanBaseline,c(2,3,4),mean,na.rm=TRUE)
#   
  
  # Transform any variables stored in log space
  resultsTimesMeanBaseline[,varsToTransform,,]<-exp(resultsTimesMeanBaseline[,varsToTransform,,])-1
  
  # Convert any variables stored in grams to kilograms
  resultsTimesMeanBaseline[,varsToTransformKG,,]<-resultsTimesMeanBaseline[,varsToTransformKG,,]/1000
  
  # Convert any plant biomass variables to total biomass
  # (not just leaf biomass)
  for (i in varsToTransformPlantBiomass){
    for (s in 1:length(sims)){
      resultsTimesMeanBaseline[s,i,,]<-resultsTimesMeanBaseline[s,i,,]/
        .GetLPRatios(resultsDir)
    }
    
  }
  
  # Sum values across the specified variables
  resultsTotalBaseline<-apply(resultsTimesMeanBaseline,c(1,3,4),sum,na.rm=TRUE)
  
  resultsTotalBaseline[resultsTotalBaseline==0]<-NA
  
  
  .Log("Processing impact results\n")
  
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
  
  # Select only the designated time steps
  allResults<-allResults[,match(times,allTimes),,,,drop=FALSE]
  
  allResults[is.na(allResults)]<-0
  
  # Calcualte the mean value across these time steps for 
  # each variable within each simulation
  resultsTimesMean<-apply(allResults,c(1,3,4,5),mean,na.rm=TRUE)
  
#   # Calcualte the mean value across simulations
#   resultsSimMean<-apply(resultsTimesMean,c(2,3,4),mean,na.rm=TRUE)
#   
  # Transform any variables stored in log space
  resultsTimesMean[,varsToTransform,,]<-exp(resultsTimesMean[,varsToTransform,,])-1
  
  # Convert any variables stored in grams to kilograms
  resultsTimesMean[,varsToTransformKG,,]<-resultsTimesMean[,varsToTransformKG,,]/1000
  
  # Convert any plant biomass variables to total biomass
  # (not just leaf biomass)
  for (i in varsToTransformPlantBiomass){
    for (s in 1:length(sims)){
      resultsTimesMean[s,i,,]<-resultsTimesMean[s,i,,]/
        .GetLPRatios(resultsDir)
    }
    
  }
  
  # Sum values across the specified variables
  resultsTotal<-apply(resultsTimesMean,c(1,3,4),sum,na.rm=TRUE)
  
  resultsTotal[resultsTotal==0]<-NA
  
  resultsPercentage <- (resultsTotal/resultsTotalBaseline)*100
  
  # Calculate average across simulations
  resultsPercentage <- apply(resultsPercentage,c(2,3),median)
  
  resultsMap<-SpatialGridDataFrame(grid = GridTopology(
    cellcentre.offset = c(min(long),min(lat)),
    cellsize = c(unique(diff(long)),unique(diff(lat))),
    cells.dim = c(length(long),length(lat))),
    data = data.frame(band1=as.vector(resultsPercentage[,ncol(resultsPercentage):1])))
  
  gridAverage <- w.median(x = resultsMap$band1,
                          w = DegreeCellAreaKM(
                            lat = coordinates(resultsMap)[,2],
                            height = resultsMap@grid@cellsize[2],
                            width = resultsMap@grid@cellsize[1]))
  
#   gridAverage <- weighted.mean(x = resultsMap$band1,
#                                w = DegreeCellAreaKM(
#                                  lat = coordinates(resultsMap)[,2],
#                                  height = resultsMap@grid@cellsize[2],
#                                  width = resultsMap@grid@cellsize[1]),na.rm=TRUE)
#   
  .Log(paste("Grid-average for ",plotName,": ",gridAverage,"\n",sep=""))
  
  .Log("Plotting\n")
  
  if(is.null(brks)){
    brks <- quantile(resultsPercentage,probs = seq(from=0,to=1,length.out=11),
                     na.rm = TRUE)
  }
  if(is.null(col)){
    cols <- brewer.pal(n = 10,name = "RdYlBu")
  } else {
    cols <- col
  }
  
  image(resultsMap,breaks = brks,col=cols,
        xaxt="n",yaxt="n",bty="n")
  
  if (map == "Africa"){
    AfricaOutline(add=TRUE)
  } else {
    plot(outline,add=TRUE)
  }
  
  invisible(dev.off())
  
  if(returnMap){
    return(resultsMap)
  }
  
}