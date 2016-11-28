
PlotPyramids <- function(resultsDir,plotName,outDir=NULL,
                         label=NULL,
                         gridSimulation=FALSE,
                        whichCells=NULL,endTimeStep=NULL,
                        numTimeSteps=12,vars=c("autotroph biomass density",
                                               "herbivore biomass density",
                                               "omnivore biomass density",
                                               "carnivore biomass density"),
                        cols=c("#1b9e77",
                                    "#66a61e",
                                    "#7570b3",
                                    "#d95f02"),
                        plotFlows=FALSE){
  
  stopifnot(length(vars)==length(cols))
  
  if ((gridSimulation) & (is.null(whichCells))){
    stop("Error, if a grid-based simulation, you must specify cells")
  }
  
  PlantBiomassVariables<-list(
    "autotroph biomass density" = TRUE,
    "carnivore density" = FALSE,
    "carnivore biomass density" = FALSE,
    "herbivore density" = FALSE,
    "herbivore biomass density" = FALSE,
    "omnivore density" = FALSE,
    "omnivore biomass density" = FALSE
  )
  
  stopifnot(all(vars %in% names(PlantBiomassVariables)))
  
  if ("SimulationControlParameters.csv" %in% dir(resultsDir)){
    initialization <- read.csv(paste(resultsDir,"/SimulationControlParameters.csv",sep=""))
  } else {
    initialization <- read.csv(paste(resultsDir,"/EcosystemModelInitialisation.csv",sep=""))
  }

  cellsize <- as.numeric(paste(initialization$Value[which(initialization$Parameter=="Grid Cell Size")]))
  
  .Log("Finding Madingley output files\n")
  if (gridSimulation){
    files<-.ListGridOutputFiles(resultsDir)
  } else {
    files<-.ListCellOuputFiles(resultsDir)
  }
  
  if (!gridSimulation){
    if(!is.null(whichCells)){
      files <- files[sapply(paste("Cell",whichCells-1,sep=""),FUN = function(x) return(grep(x,files)))]
    }
    
    # Find the unique cells in these simulations
    cells.re<-regexpr("Cell[0-9]+",files)
    cells<-as.list(unique(substr(files,cells.re,cells.re+
                                   attr(cells.re,"match.length")-1)))
    
  }
  
  if (gridSimulation) {
    # Find the simulation numbers
    sims.re<-regexpr("_[0-9]+",files)
    sims<-as.list(unique(substr(files,sims.re,sims.re+
                                  attr(sims.re,"match.length")-1)))
    
  } else {
    # Find the simulation numbers
    sims.re<-regexpr("_[0-9]+_",files)
    sims<-as.list(unique(substr(files,sims.re,sims.re+
                                  attr(sims.re,"match.length")-1)))
  }
  
  if(is.null(label)){
    label<-unique(substr(files,1,sims.re-1))
    print(label)
    stopifnot(length(label)==1)
    label<-label[1]
  } else {
    label <- paste("BasicOutputs_",label,sep="")
  }
  
  if (!gridSimulation) .Log(paste("Found results for ",length(cells)," cells\n",sep=""))
  .Log(paste("Found results for ",length(sims)," simulations\n",sep=""))
  
  .Log("Getting basic information about simulations\n")
  if (gridSimulation){
    sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],
                    ".nc",sep="")
  } else {
    sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],cells[1],
                    ".nc",sep="")
  }
  
  data<-open.sds(sds.path)
  allTimes<-get.sds(data,"Time step")
  if (is.null(endTimeStep)) endTimeStep <- tail(allTimes,1)
  times <- (endTimeStep-numTimeSteps+1):endTimeStep
  
  .Log("Initializing plot\n")
  if (gridSimulation){
    dims <- .plotDims(length(whichCells))
  } else {
    dims<-.plotDims(length(cells))
  }
  
  if(!is.null(outDir)){
    pdf(paste(outDir,plotName,".pdf",sep=""),
        width = dims$width,height = dims$height)
  }
  
  if (gridSimulation){
    par(mfrow=.gridArrange(length(whichCells)))
  } else {
    par(mfrow=.gridArrange(length(cells)))
  }
  
  .Log("Plotting\n")
  lapply(cells,FUN=function(cell){
    
    # Create a list of matrices to hold the results for each specified variable
    allResults<-list()
    for (i in 1:length(vars)){
      allResults[[i]]<-matrix(data = NA,nrow = length(sims),
                              ncol = length(allTimes))
    }
    names(allResults)<-vars
    
    # Loop over simulations in the ensemble
    s<-1
    for (sim in sims){
      if (gridSimulation){
        sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,sep="")
      } else {
        sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,cell,
                        ".nc",sep="")
      }
      
      data<-open.sds(sds.path)
      
      # Populate the results matrices
      for (var in vars){
        test <- get.sds(data,var)
        print("Input data dimensions:")
        print(dim(test))
        print("Destination dimensions:")
        print(dim(allResults[var][[1]][s,]))
        stop()
        allResults[var][[1]][s,]<-get.sds(data,var)
      }
      s<-s+1
    }
    
    for (var in vars){
      if (PlantBiomassVariables[var][[1]]){
        lp.ratio <- .GetLPRatios(resultsDir)[as.integer(gsub("Cell","",cell))+1]
        allResults[[var]] <- allResults[[var]]/lp.ratio
      }
    }
    
    resultsTimesMean <- lapply(allResults,function(x){
      return(apply(x[,times,drop=FALSE],MARGIN=1,FUN=mean,na.rm=TRUE))
    })
    
    resultsSimMean <- lapply(resultsTimesMean,FUN=mean,na.rm=TRUE)
    
    y.breaks <- seq(from=0,to=1,length.out=length(vars)+1)
    
    temp.df<-data.frame(y1=y.breaks[1:length(y.breaks)-1],y2=y.breaks[2:length(y.breaks)])
    
    yVals <- as.list(as.data.frame(t(temp.df)))
    
    maxBiomass <- max(unlist(resultsSimMean))
    relWidths <- lapply(resultsSimMean,function(x){
      return(log10(x)/log10(maxBiomass))
    })
    
    par(mar=c(1,1,1,8.5))
    plot.new()
    
    mapply(FUN = function(width,y,col,biomass){
      x1 <- (1-width)/2
      x2 <- 1 - x1
      rect(xleft = x1,ybottom = y[1],xright = x2,ytop = y[2],col=col,border=NA)
      text(0.5,mean(y),paste(round(biomass/1000,1),"t"))
    },width=relWidths,y=yVals,col=as.list(cols),biomass=resultsSimMean)
    
    if (plotFlows){
      
      flows <- list()
      s <- 1
      for (sim in sims){
        path <- paste(gsub("BasicOutputs","TrophicFlows",label),sim,cell,".txt",sep="")
        
        flows[[s]] <- read.table(paste(resultsDir,"/",path,sep=""),header=TRUE)
        
        flows[[s]] <- flows[[s]][(flows[[s]]$time_step %in% times),]
        
        s <- s+1
      }
      
      allSims <- Reduce(function(...) return(merge(...,all=TRUE)),flows)
      
      allSims <- aggregate(mass_eaten_g ~ Longitude + Latitude + time_step + fromIndex + toIndex,
                           FUN=mean,data=allSims)
      
      allSims <- aggregate(mass_eaten_g ~ Longitude + Latitude + fromIndex + toIndex,
                           FUN=sum,data=allSims)
      
      if ((length(unique(allSims$Longitude)) != 1) || 
          (length(unique(allSims$Latitude)) != 1))
      stop("Error: this output file is for more than one cell")
      
      allSims$mass_eaten_g<-allSims$mass_eaten_g/DegreeCellAreaKM(
        lat = as.numeric(paste(allSims$Latitude))+cellsize/2,
        height = cellsize,width = cellsize)
      
      allSims$mass_eaten_g[which(allSims$toIndex==0)]<-allSims$mass_eaten_g[which(allSims$toIndex==0)]/
        .GetLPRatios(resultsDir)[as.integer(gsub("Cell","",cell))+1]
      
      allSims$proportion <- allSims$mass_eaten_g/1000/
        unlist(resultsSimMean[as.numeric(paste(allSims$fromIndex))+1])
      
    }
    
    legend(x = 1.1,y = 0.5,
           legend = gsub(" biomass density","",vars),
           fill = cols,xpd=TRUE,bty = "n")
    
#     
#     yearAvgs<-lapply(allResults,FUN = function(resultsMatrix){
#       return(
#         apply(resultsMatrix,1,FUN = function(simResults){
#           tapply(simResults,y,mean)
#         })
#       )
#     })
#     
#     # For each variable get the median and confidence limits across simulations
#     summaryValues<-lapply(yearAvgs,FUN = function(x){
#       medianValues<-apply(x,1,median)
#       upperValues<-apply(x,1,quantile,
#                          probs=c(1-(1-(confidenceInterval/100))/2),
#                          na.rm=TRUE)
#       lowerValues<-apply(x,1,quantile,
#                          probs=c((1-(confidenceInterval/100))/2),
#                          na.rm=TRUE)
#       
#       medianValues[medianValues==0]<-NA
#       upperValues[upperValues==0]<-NA
#       upperValues[is.na(medianValues)]<-NA
#       lowerValues[lowerValues==0]<-NA
#       lowerValues[is.na(medianValues)]<-NA
#       
#       return(list(median=medianValues,
#                   lower=lowerValues,
#                   upper=upperValues))
#     })
#     
#     maxVal<-max(unlist(lapply(summaryValues,function(x){
#       return(max(x$upper,na.rm=TRUE))}
#     )))
#     minVal<-min(unlist(lapply(summaryValues,function(x){
#       return(min(x$lower,na.rm=TRUE))
#     })))
#     
#     par(mgp=c(1.6,0.3,0))
#     plot(999999999999999,999999999999999,xlim=range(years),ylim=c(minVal,maxVal),log="y",
#          xlab=xlab,yaxt="n",ylab=NA)
#     par(mgp=c(2.6,0.3,0))
#     axis(2)
#     title(ylab=ylab)
#     
#     if (plotConfidence){
#       mapply(FUN=function(summaryValues,cols){
#         X.Vec<-c(years,max(years),rev(years),min(years))
#         Y.Vec<-c(summaryValues$lower,tail(summaryValues$upper,1),
#                  rev(summaryValues$upper),summaryValues$lower[1])
#         X.Vec<-X.Vec[!is.na(Y.Vec)]
#         Y.Vec<-Y.Vec[!is.na(Y.Vec)]
#         polygon(X.Vec,Y.Vec,col=paste(cols,"33",sep=""),border=NA)
#       },summaryValues,cols)
#     }
#     
#     mapply(FUN = function(summaryValues,cols){
#       points(years,summaryValues$median,type="l",col=cols)
#     },summaryValues,as.list(cols))
#     
  })
  
  if(!is.null(outDir)) invisible(dev.off())
}