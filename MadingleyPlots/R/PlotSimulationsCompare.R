
PlotSimulationsCompare <- function(resultsDirBaseline,resultsDir,plotName,outDir,
                               labelBaseline=NULL,labels=NULL,xVals=NULL,
                               vars=c("autotroph biomass density",
                                      "herbivore biomass density",
                                      "omnivore biomass density",
                                      "carnivore biomass density"),
                               cols=c("#1b9e77",
                                      "#66a61e",
                                      "#7570b3",
                                      "#d95f02"),
                               ylab=plotName,xlab="Impact level",
                               confidenceInterval=95,
                               endTimeStepBaseline=NULL,numTimeStepsBaseline=12,
                               endTimeStep=NULL,numTimeSteps=12,
                               whichCells=NULL,
                               returnResults = FALSE,
                               pairedSimulations=FALSE){
  
  PlantBiomassVariables<-list(
    "autotroph biomass density" = TRUE,
    "carnivore density" = FALSE,
    "carnivore biomass density" = FALSE,
    "herbivore density" = FALSE,
    "herbivore biomass density" = FALSE,
    "omnivore density" = FALSE,
    "omnivore biomass density" = FALSE,
    "Generalist herbivoreBiomass" = FALSE,
    "Primary herbivoreBiomass" = FALSE,
    "Secondary herbivoreBiomass" = FALSE,
    "Plantation herbivoreBiomass" = FALSE,
    "Cropland herbivoreBiomass" = FALSE,
    "Pasture herbivoreBiomass" = FALSE,
    "Urban herbivoreBiomass" = FALSE,
    "Generalist omnivoreBiomass" = FALSE,
    "Primary omnivoreBiomass" = FALSE,
    "Secondary omnivoreBiomass" = FALSE,
    "Plantation omnivoreBiomass" = FALSE,
    "Cropland omnivoreBiomass" = FALSE,
    "Pasture omnivoreBiomass" = FALSE,
    "Urban omnivoreBiomass" = FALSE,
    "Generalist carnivoreBiomass" = FALSE,
    "Primary carnivoreBiomass" = FALSE,
    "Secondary carnivoreBiomass" = FALSE,
    "Plantation carnivoreBiomass" = FALSE,
    "Cropland carnivoreBiomass" = FALSE,
    "Pasture carnivoreBiomass" = FALSE,
    "Urban carnivoreBiomass" = FALSE,
    "Generalist herbivoreDensity" = FALSE,
    "Primary herbivoreDensity" = FALSE,
    "Secondary herbivoreDensity" = FALSE,
    "Plantation herbivoreDensity" = FALSE,
    "Cropland herbivoreDensity" = FALSE,
    "Pasture herbivoreDensity" = FALSE,
    "Urban herbivoreDensity" = FALSE,
    "Generalist omnivoreDensity" = FALSE,
    "Primary omnivoreDensity" = FALSE,
    "Secondary omnivoreDensity" = FALSE,
    "Plantation omnivoreDensity" = FALSE,
    "Cropland omnivoreDensity" = FALSE,
    "Pasture omnivoreDensity" = FALSE,
    "Urban omnivoreDensity" = FALSE,
    "Generalist carnivoreDensity" = FALSE,
    "Primary carnivoreDensity" = FALSE,
    "Secondary carnivoreDensity" = FALSE,
    "Plantation carnivoreDensity" = FALSE,
    "Cropland carnivoreDensity" = FALSE,
    "Pasture carnivoreDensity" = FALSE,
    "Urban carnivoreDensity" = FALSE,
    "Primary Landuse Affinity" = FALSE,
    "Secondary Landuse Affinity" = FALSE,
    "Plantation Landuse Affinity" = FALSE,
    "Cropland Landuse Affinity" = FALSE,
    "Pasture Landuse Affinity" = FALSE,
    "Urban Landuse Affinity" = FALSE
    
  )
  
  stopifnot(all(vars %in% names(PlantBiomassVariables)))
  
  .Log("Finding Madingley output files\n")
  files<-.ListCellOuputFiles(resultsDir)
  if (!is.null(labels)){
    files <- files[as.logical(apply(do.call("cbind",lapply(
      as.list(labels),function(x) return(grepl(x,files)))),
      1,max))]
  }
  
  if(!is.null(whichCells)){
    files <- files[which(as.logical(apply(do.call('cbind',lapply(
      as.list(whichCells),function(x) return(grepl(paste(
        "Cell",x-1,sep=""),files)))),1,max)))]
  }
  
  filesBaseline <- .ListCellOuputFiles(resultsDirBaseline)
  if (!is.null(labelBaseline)){
    filesBaseline <- filesBaseline[grepl(labelBaseline,filesBaseline)]
  }

  if(!is.null(whichCells)){
    filesBaseline <- filesBaseline[which(as.logical(apply(do.call('cbind',lapply(
      as.list(whichCells),function(x) return(grepl(paste(
        "Cell",x-1,sep=""),filesBaseline)))),1,max)))]
  }
  
  # Find the simulation numbers
  sims.reBaseline<-regexpr("_[0-9]+_",filesBaseline)
  simsBaseline<-as.list(unique(substr(filesBaseline,sims.reBaseline,sims.reBaseline+
                                        attr(sims.reBaseline,"match.length")-1)))
  .Log(paste("Found baseline results for ",length(simsBaseline)," simulations\n",sep=""))
  
  sims.re<-regexpr("_[0-9]+_",files)
  sims<-as.list(unique(substr(files,sims.re,sims.re+
                                attr(sims.re,"match.length")-1)))
  .Log(paste("Found impact results for ",length(sims)," simulations\n",sep=""))
  
  if (length(sims) != length(simsBaseline)){
    if (pairedSimulations){
      stop("Error: number of simulations under impact not equal to number in baseline")
    } else {
      .Log("Warning: number of simulations under impact not equal to number in baseline\n")
    }
    
  }
  
  # Find the unique cells in these simulations
  cells.re<-regexpr("Cell[0-9]+",files)
  cells<-as.list(unique(substr(files,cells.re,cells.re+
                                 attr(cells.re,"match.length")-1)))
  .Log(paste("Found results for ",length(cells)," cells\n",sep=""))
  
  if(is.null(labelBaseline)){
    labelBaseline<-unique(substr(filesBaseline,1,sims.reBaseline-1))
    stopifnot(length(labelBaseline)==1)
    labelBaseline<-labelBaseline[1]
  } else {
    labelBaseline <- paste("BasicOutputs_",labelBaseline,sep="")
  }
  
  if(is.null(labels)){
    labels<-unique(substr(files,1,sims.re-1))
    stopifnot(length(labels)==1)
    labels<-labels[1]
  } else {
    labels <- paste("BasicOutputs_",labels,sep="")
  }
  
  if(is.null(xVals)){
    xVals <- 1:length(labels)
  }
  
  .Log("Getting basic information about simulations\n")
  sds.path<-paste("msds:nc?file=",resultsDir,"/",labelBaseline,simsBaseline[1],cells[1],
                  ".nc",sep="")
  data<-open.sds(sds.path)
  allTimes<-get.sds(data,"Time step")
  if (is.null(endTimeStep)) endTimeStep <- tail(allTimes,1)
  times <- (endTimeStep-numTimeSteps+1):endTimeStep
  
  .Log("Initializing plot\n")
  dims<-.plotDims(length(cells))
  pdf(paste(outDir,plotName,".pdf",sep=""),
      width = dims$width,height = dims$height)
  
  par(mfrow=.gridArrange(length(cells)))
  par(mar=c(2.5,2.5,0.5,0.5))
  par(tck=-0.01)
  par(mgp=c(1.6,0.2,0))
  par(las=1)
  
  .Log("Plotting\n")
  ret <- list()
  lapply(cells,FUN=function(cell){
    
    cat(paste("Processing Cell: ",cell,"\n",sep=""))
    
    # Create matrices to hold the results for each specified variable
    allResults<-list()
    for (i in 1:length(vars)){
      allResults[[i]]<-data.frame(matrix(data = NA,nrow = length(sims),
                                         ncol = length(labels)))
      names(allResults[[i]]) <- xVals
    }
    names(allResults)<-vars
    
    if (pairedSimulations){
      
      # Loop over simulations in the ensemble
      s<-1
      for (sim in sims){
        
        for (var in vars){
          if (PlantBiomassVariables[var][[1]]){
            lp.ratio <- .GetLPRatios(resultsDir)[as.integer(gsub("Cell","",cell))+1]
          }
          
          sds.path.bl<-paste("msds:nc?file=",resultsDir,"/",labelBaseline,sim,cell,
                             ".nc",sep="")
          data.bl<-open.sds(sds.path.bl)
          timeSeries.bl<-get.sds(data.bl,var)[times]
          
          if (PlantBiomassVariables[var][[1]]){
            timeSeries.bl <- timeSeries.bl/lp.ratio
          }
          
          timeAvg.bl <- mean(timeSeries.bl)
          
          for (l in 1:length(labels)){
            sds.path<-paste("msds:nc?file=",resultsDir,"/",labels[l],sim,cell,".nc",sep="")
            data<-open.sds(sds.path)
            timeSeries<-get.sds(data,var)[times]
            
            if (PlantBiomassVariables[var][[1]]){
              timeSeries <- timeSeries/lp.ratio
            }
            
            timeAvg <- mean(timeSeries)
            
            allResults[[var]][s,l] <- (timeAvg/timeAvg.bl)*100
            
          }
          
          
          
        }
        
        s <- s+1
        
        
      }
      
    } else {
      stop("Independent baseline and impact simulations not supported at present")
    }
    
    if (length(vars)==1){
      offsets <- 0
    } else {
      offsets <- seq(from=-0.25,to=0.25,length.out=length(vars))
    }
    
    combinedResults <- do.call('cbind',allResults)
    lowers <- apply(combinedResults,2,quantile,probs=((1-(confidenceInterval/100))/2),na.rm=TRUE)
    uppers <- apply(combinedResults,2,quantile,probs=1-((1-(confidenceInterval/100))/2),na.rm=TRUE)
    yMax <- max(uppers,na.rm=TRUE)
    yMin <- min(lowers,na.rm=TRUE)
    
    errbar(-9e99,-9e99,-9e99,-9e99,xlim=c(1+min(offsets),
                                          dim(allResults[[1]])[2]+max(offsets)),
           ylim=c(yMin,yMax),xaxt="n",xlab=NA,ylab=ylab)
    par(mgp=c(1,0.2,0))
    axis(1,at=1:dim(allResults[[1]])[2],labels=xVals)
    title(xlab=xlab)
    par(mgp=c(1.6,0.2,0))
    v <- 1
    for (var in vars){
      errbar((1:dim(allResults[[1]])[2])+offsets[v],
             apply(allResults[[var]],2,median,na.rm=TRUE),
             apply(allResults[[var]],2,quantile,probs=((1-(confidenceInterval/100))/2),na.rm=TRUE),
             apply(allResults[[var]],2,quantile,probs=1-((1-(confidenceInterval/100))/2),na.rm=TRUE),
             add=TRUE,col=cols[v],errbar.col=cols[v])
      
      v <- v+1
      
    }
    abline(h=100,lty=2,col="#00000033")
    
  })
  
  
  invisible(dev.off())
  
}
