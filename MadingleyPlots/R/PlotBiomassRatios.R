
PlotBiomassRatios <- function(resultsDir,plotName="BiomassRatios",
                              label=NULL,
                              outDir=NULL,cellLong,cellLat,
                              endTimeStep=NULL,
                              numTimeSteps=12,
                              vars=c("herbivorebiomass density",
                                     "omnivorebiomass density",
                                     "carnivorebiomass density"),
                              labels=NULL,
                              cols=c("#66a61e",
                                     "#7570b3",
                                     "#d95f02"),
                              ylims=NULL,
                              returnResults=FALSE){
  
  stopifnot(length(vars)==length(cols))
  
  LogVariables<-list(
    "Biomass density" = TRUE,
    "autotrophbiomass density" = TRUE,
    "carnivorebiomass density" = TRUE,
    "herbivorebiomass density" = TRUE,
    "omnivorebiomass density" = TRUE
    
  )
  
  GramsBiomassVariables<-list(
    "Biomass density" = TRUE,
    "autotrophbiomass density" = TRUE,
    "carnivorebiomass density" = TRUE,
    "herbivorebiomass density" = TRUE,
    "omnivorebiomass density" = TRUE
    
  )
  
  PlantBiomassVariables<-list(
    "Biomass density" = FALSE,
    "autotrophbiomass density" = TRUE,
    "carnivorebiomass density" = FALSE,
    "herbivorebiomass density" = FALSE,
    "omnivorebiomass density" = FALSE
    
  )
  
  stopifnot(all(vars %in% names(LogVariables)))
  stopifnot(all(vars %in% names(GramsBiomassVariables)))
  stopifnot(all(vars %in% names(PlantBiomassVariables)))
  
  varsToTransform<-which(unlist(LogVariables[match(c("autotrophbiomass density",vars),names(LogVariables))]))
  varsToTransformKG<-which(unlist(GramsBiomassVariables[
    match(c("autotrophbiomass density",vars),names(GramsBiomassVariables))]))
  varsToTransformPlantBiomass<-which(unlist(PlantBiomassVariables[match(
    c("autotrophbiomass density",vars),names(PlantBiomassVariables))]))
  
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
  
  stopifnot((min(long) <= cellLong) & (cellLong <= max(long)))
  stopifnot((min(lat) <= cellLat) & (cellLat <= max(lat)))
  
  whichLongCell <- order(abs(cellLong - long))[1]
  whichLatCell <- order(abs(cellLat - lat))[1]

  allTimes <- get.sds(data,"Time step")
  if (is.null(endTimeStep)) endTimeStep <- tail(allTimes,1)
  times <- (endTimeStep-numTimeSteps+1):endTimeStep
  
  .Log("Initializing plot\n")
  if (!is.null(outDir)){
    pdf(paste(outDir,plotName,".pdf",sep=""),
        width = 8.5/2.54,height = 8.5/2.54)
  }
  
  par(las=1)
  par(tck=-0.01)
  par(mgp=c(1.2,0.2,0))
  par(mar=c(2.6,2.4,0.2,0.2))
  
  allResults <- array(data = NA,dim = c(length(vars)+1,length(sims),length(allTimes),length(long),length(lat)))
  
  # Loop over simulations in the ensemble
  s <- 1
  for (sim in sims){
    sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,sep="")
    
    data<-open.sds(sds.path)
    
    v <- 1
    for (var in c("autotrophbiomass density",vars)){
      allResults[v,s,,,]<- get.sds(data,var)
      v <- v+1
    }
    
    s <- s+1
  }
  
  # Select only the designated time steps
  allResults<-allResults[,,match(times,allTimes),,,drop=FALSE]
  
  allResults[is.na(allResults)]<-0
  
  # Calcualte the mean value across the selected time steps
  resultsTimesMean <- apply(allResults,MARGIN = c(1,2,4,5),FUN = mean,na.rm=TRUE)
  
  # Transform any variables stored in log space
  resultsTimesMean[varsToTransform,,,]<-exp(resultsTimesMean[varsToTransform,,,])-1
  
  # Convert any variables stored in grams to kilograms
  resultsTimesMean[varsToTransformKG,,,]<-resultsTimesMean[varsToTransformKG,,,]/1000
  
  # Convert any plant biomass variables to total biomass
  # (not just leaf biomass)
  for (i in varsToTransformPlantBiomass){
    for (s in 1:length(sims)){
      resultsTimesMean[i,s,,]<-resultsTimesMean[i,s,,]/
        .GetLPRatios(resultsDir)
    }
    
  }
  
  # Calculate log10 ratios of variables to autotroph biomass
  for (v in 2:dim(resultsTimesMean)[1]){
    resultsTimesMean[v,,,] <- log10(resultsTimesMean[v,,,]/resultsTimesMean[1,,,])
  }
  
  # Calculate medians and confidence limits of biomass ratios across simulations
  medians <- apply(resultsTimesMean,MARGIN = c(1,3,4),FUN = median)
  uppers <- apply(resultsTimesMean,MARGIN = c(1,3,4),FUN = quantile,probs=0.975,na.rm=TRUE)
  lowers <- apply(resultsTimesMean,MARGIN = c(1,3,4),FUN = quantile,probs=0.025,na.rm=TRUE)
  
  se <- function(x) sqrt(var(na.omit(x))/length(na.omit(x)))
  
  means <- apply(resultsTimesMean, MARGIN = c(1,3,4), FUN = mean)
  ses <- apply(resultsTimesMean, MARGIN = c(1,3,4), FUN = se)
  
  # Select the specified cell
  medians <- medians[,whichLongCell,whichLatCell][2:(length(vars)+1)]
  uppers <- uppers[,whichLongCell,whichLatCell][2:(length(vars)+1)]
  lowers <- lowers[,whichLongCell,whichLatCell][2:(length(vars)+1)]
  
  means <- means[,whichLongCell,whichLatCell][2:(length(vars)+1)]
  ses <- ses[,whichLongCell,whichLatCell][2:(length(vars)+1)]
  
  if(is.null(labels)){
    labels <- vars
  }
  
  names(means) <- labels
  names(ses) <- labels
  
  errbar(1:3,medians,uppers,lowers,ylim=ylims,xlim=c(0.5,length(vars)+0.5),
         pch=23,col=cols,errbar.col=cols,cex=2,xaxt="n",xlab=NA,
         ylab=expression(paste("Log"[10]," Ratio to Autotroph Biomass",sep="")))
  points(1:3,medians,pch=23,col=cols,bg="white",cex=2)
  axis(1,at=1:length(vars),labels=labels)
  
  if (!is.null(outDir)){
    invisible(dev.off())
  }
  
  if (returnResults) {
    return(list(mean=means,se=ses))
  }
  
}