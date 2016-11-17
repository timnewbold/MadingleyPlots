
PlotMassDensity <- function(resultsDir,plotName = "MassDensity",
                            outDir=NULL,
                            label=NULL,
                            whichCells=NULL,endTimeStep=NULL,
                            numTimeSteps=12,
                            vars=c("herbivore abundance",
                                   "omnivore abundance",
                                   "carnivore abundance"),
                            cols=c("#66a61e",
                                   "#7570b3",
                                   "#d95f02"),
                            xlims = NULL,
                            returnResults=FALSE){
  
  initialization <- read.csv(paste(resultsDir,"/SimulationControlParameters.csv",sep=""))
  cellsize <- as.numeric(paste(initialization$Value[which(initialization$Parameter=="Grid Cell Size")]))
  
  locations <- read.csv(paste(resultsDir,"/SpecificLocations.csv",sep=""))
  
  cohortDefs <- read.csv(paste(resultsDir,"/CohortFunctionalGroupDefinitions.csv",sep=""))
  maxPossibleMass <- max(cohortDefs$PROPERTY_Maximum.mass)
  minPossibleMass <- min(cohortDefs$PROPERTY_Minimum.mass)
  
  .Log("Finding Madingley mass-bins output files\n")
  files <- .ListMassBinsFiles(resultsDir)
  
  if(!is.null(whichCells)){
    files <- files[sapply(paste("Cell",whichCells-1,sep=""),FUN = function(x) return(grep(x,files)))]
  }
  
  # Find the unique cells in these simulations
  cells.re<-regexpr("Cell[0-9]+",files)
  cells<-as.list(unique(substr(files,cells.re,cells.re+
                                 attr(cells.re,"match.length")-1)))
  
  # Find the simulation numbers
  sims.re<-regexpr("_[0-9]+_",files)
  sims<-as.list(unique(substr(files,sims.re,sims.re+
                                attr(sims.re,"match.length")-1)))
  
  if(is.null(label)){
    label<-unique(substr(files,1,sims.re-1))
    print(label)
    stopifnot(length(label)==1)
    label<-label[1]
  } else {
    label <- paste("MassBinsOutputs_",label,sep="")
  }
  
  .Log(paste("Found results for ",length(cells)," cells\n",sep=""))
  .Log(paste("Found results for ",length(sims)," simulations\n",sep=""))
  
  .Log("Getting basic information about simulations\n")
  sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],cells[1],
                  ".nc",sep="")
  data<-open.sds(sds.path)
  allTimes<-get.sds(data,"Time step")
  if (is.null(endTimeStep)) endTimeStep <- tail(allTimes,1)
  times <- (endTimeStep-numTimeSteps+1):endTimeStep
  
  massBins <- get.sds(data,"Mass bin")
  
  massBinsMidPoints <- 10^(log10(c(minPossibleMass,massBins[-1]))+diff(log10(c(
    minPossibleMass,massBins[-1],maxPossibleMass)))/2)
  
  latitudes <- locations$Latitude
  if(!is.null(whichCells)){
    latitudes <- latitudes[whichCells]
  }
  cell_areas <- DegreeCellAreaKM(lat = latitudes,height = cellsize,width = cellsize)
  
  names(cell_areas) <- cells
  
  if (is.null(xlims)){
    xlims <- range(massBinsMidPoints)
  }
  
  .Log("Initializing plot\n")
  if(!is.null(outDir)){
    pdf(paste(outDir,plotName,".pdf",sep=""),
        width = 17.5/2.54,height = (5/2.54)*length(cells))
  }
  
  par(mfrow=c(length(cells),3))
  par(las=1)
  par(tck=-0.01)
  par(mar=c(2.8,3.3,0.2,0.2))
  
  .Log("Plotting\n")
  
  ret <- lapply(cells,FUN=function(cell){
    
    # Create a list of matrices to hold the results for each specified variable
    allResults<-list()
    for (i in 1:length(vars)){
      allResults[[i]]<-array(data = NA,dim = c(length(sims),length(massBins),length(allTimes)))
    }
    names(allResults)<-vars
    
    # Loop over simulations in the ensemble
    s<-1
    for (sim in sims){
      sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,cell,
                      ".nc",sep="")
      
      data<-open.sds(sds.path)
      
      # Populate the results matrices
      for (var in vars){
        allResults[var][[1]][s,,]<-exp(get.sds(data,paste(
          "Log ",var," in mass bins",sep="")))/cell_areas[cell]
      }
      s<-s+1
    }
    
    resultsTimesMean <- lapply(allResults,function(x){
      return(apply(x[,,times,drop=FALSE],MARGIN=c(1,2),FUN=mean,na.rm=TRUE))
    })
    
    resultsSimMean <- lapply(resultsTimesMean,function(x){
      return(apply(x,MARGIN=2,FUN=mean,na.rm=TRUE))
    })
    
    v <- 1
    r <- list()
    for (var in vars){
      par(mgp=c(1.2,0.2,0))
      plot(massBinsMidPoints,resultsSimMean[[var]],log="xy",type="l",
           col=cols[v],xlim=xlims,xlab="Current body mass (g)",
           yaxt="n",ylab=NA)
      par(mgp=c(2.5,0.2,0))
      axis(2)
      title(ylab=Hmisc::capitalize(var))
      v <- v+1
      r[[var]] <- data.frame(mass=massBinsMidPoints,
                             density=resultsSimMean[[var]])
    }
    
    return(r)
    
  })
  
  if(!is.null(outDir)) invisible(dev.off())
  
  if (returnResults){
    return(ret)
  }
  
}