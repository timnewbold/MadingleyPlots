
PlotTemporal<-function(resultsDir,plotName,outDir=NULL,
                       label=NULL,
                       gridSimulation=FALSE,
                       whichCells = NULL,
                       vars=c("autotroph biomass density",
                              "herbivore biomass density",
                              "omnivore biomass density",
                              "carnivore biomass density"),
                       cols=c("#1b9e77",
                              "#66a61e",
                              "#7570b3",
                              "#d95f02"),
                       xlab="Years",ylab=plotName,
                       plotConfidence=TRUE,
                       confidenceInterval=95,
                       returnResults = FALSE){
  
  stopifnot(length(vars)==length(cols))
  
  if ((gridSimulation) & (is.null(whichCells))){
    stop("Error, if a grid-based simulation, you must specify cells")
  }
  
  if (gridSimulation){
    vars <- gsub(" biomass density","biomass density",vars)
  }
  
  PlantBiomassVariables<-list(
    "autotroph biomass density" = TRUE,
    "carnivore density" = FALSE,
    "carnivore biomass density" = FALSE,
    "herbivore density" = FALSE,
    "herbivore biomass density" = FALSE,
    "omnivore density" = FALSE,
    "omnivore biomass density" = FALSE,
    "Mean Trophic Level" = FALSE,
    "Max Trophic Index" = FALSE,
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
  
  if (gridSimulation){
    names(PlantBiomassVariables) <- lapply(
      names(PlantBiomassVariables),function(x) return(
        gsub(" biomass density","biomass density",x)))
  }
  
  if (gridSimulation){
    
    LogVariables<-list(
      "Abundance density" = TRUE,
      "Biomass density" = TRUE,
      "autotrophbiomass density" = TRUE,
      "carnivoreabundance density" = TRUE,
      "carnivorebiomass density" = TRUE,
      "herbivoreabundance density" = TRUE,
      "herbivorebiomass density" = TRUE,
      "omnivoreabundance density" = TRUE,
      "omnivorebiomass density" = TRUE,
      "Mean Trophic Level" = FALSE,
      "Max Trophic Index" = FALSE
      
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
      "omnivorebiomass density" = TRUE,
      "Mean Trophic Level" = FALSE,
      "Max Trophic Index" = FALSE
    )
    
    stopifnot(all(vars %in% names(LogVariables)))
    stopifnot(all(vars %in% names(GramsBiomassVariables)))
  }
  
  stopifnot(all(vars %in% names(PlantBiomassVariables)))
  
  if ("SimulationControlParameters.csv" %in% dir(resultsDir)){
    initialization <- read.csv(paste(resultsDir,"/SimulationControlParameters.csv",sep=""))
  } else {
    initialization <- read.csv(paste(resultsDir,"/EcosystemModelInitialisation.csv",sep=""))
  }
  
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
  times<-get.sds(data,"Time step")
  y<-times/12
  y<-ceiling(y)
  years<-unique(y)
  
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
  
  if (gridSimulation) cells <- as.list(1:length(whichCells))
  
  par(mfrow=.gridArrange(length(cells)))
  par(mar=c(3,3.5,0.5,8.5))
  par(tck=-0.01)
  par(las=1)
  
  .Log("Plotting\n")
  ret <- list()
  lapply(cells,FUN=function(cell){
    
    # Create matrices to hold the results for each specified variable
    allResults<-list()
    for (i in 1:length(vars)){
      allResults[[i]]<-matrix(data = NA,nrow = length(sims),
                              ncol = length(times))
    }
    names(allResults)<-vars
    
    # Loop over simulations in the ensemble
    s<-1
    for (sim in sims){
      if (gridSimulation){
        sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,".nc",sep="")
      } else {
        sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,cell,
                        ".nc",sep="")
      }
      
      data<-open.sds(sds.path)
      
      # Populate the results matrices
      for (var in vars){
        if (gridSimulation){
          allResults[var][[1]][s,]<-get.sds(data,var)[,whichCells[[cell]][1],whichCells[[cell]][2]]
        } else {
          allResults[var][[1]][s,]<-get.sds(data,var)
        }
      }
      s<-s+1
    }
    
    for (var in vars){
      if (gridSimulation){
        if (LogVariables[var][[1]]){
          allResults[[var]] <- exp(allResults[[var]])-1
        }
        if (GramsBiomassVariables[var][[1]]){
          allResults[[var]] <- allResults[[var]]/1000.0
        }
      }
      if (PlantBiomassVariables[var][[1]]){
        lp.ratio <- .GetLPRatios(resultsDir)[as.integer(gsub("Cell","",cell))+1]
        allResults[[var]] <- allResults[[var]]/lp.ratio
      }
    }
    
    yearAvgs<-lapply(allResults,FUN = function(resultsMatrix){
      return(
        apply(resultsMatrix,1,FUN = function(simResults){
          tapply(simResults,y,mean)
        })
      )
    })
    
    # For each variable get the median and confidence limits across simulations
    summaryValues<-lapply(yearAvgs,FUN = function(x,ret){
      medianValues<-apply(x,1,median)
      upperValues<-apply(x,1,quantile,
                         probs=c(1-(1-(confidenceInterval/100))/2),
                         na.rm=TRUE)
      lowerValues<-apply(x,1,quantile,
                         probs=c((1-(confidenceInterval/100))/2),
                         na.rm=TRUE)
      
      medianValues[medianValues==0]<-NA
      upperValues[upperValues==0]<-NA
      upperValues[is.na(medianValues)]<-NA
      lowerValues[lowerValues==0]<-NA
      lowerValues[is.na(medianValues)]<-NA
      
      return(list(median=medianValues,
                  lower=lowerValues,
                  upper=upperValues))
    })
    
    maxVal<-max(unlist(lapply(summaryValues,function(x){
      return(max(x$upper,na.rm=TRUE))}
      )))
    minVal<-min(unlist(lapply(summaryValues,function(x){
      return(min(x$lower,na.rm=TRUE))
    })))
    
    par(mgp=c(1.6,0.3,0))
    plot(999999999999999,999999999999999,xlim=range(years),
         ylim=c(minVal,maxVal),log="y",
         xlab=xlab,yaxt="n",ylab=NA,bty="l")
    par(mgp=c(2.6,0.3,0))
    axis(2)
    title(ylab=ylab)
    
    if (plotConfidence){
      mapply(FUN=function(summaryValues,cols){
        X.Vec<-c(years,max(years),rev(years),min(years))
        Y.Vec<-c(summaryValues$lower,tail(summaryValues$upper,1),
                 rev(summaryValues$upper),summaryValues$lower[1])
        X.Vec<-X.Vec[!is.na(Y.Vec)]
        Y.Vec<-Y.Vec[!is.na(Y.Vec)]
        polygon(X.Vec,Y.Vec,col=paste(cols,"33",sep=""),border=NA)
      },summaryValues,cols)
    }
    
    mapply(FUN = function(summaryValues,cols){
      points(years,summaryValues$median,type="l",col=cols)
    },summaryValues,as.list(cols))
    
    if (returnResults){
      ret[[cell]] <<- summaryValues
    }
    
    legend(x = range(years)[2]+diff(range(years))*0.07,
           y = 10^(log10(minVal)+0.7*(log10(maxVal)-log10(minVal))),
           legend = gsub(" biomass density","",vars),
           lty=1,col = cols,xpd=TRUE,bty = "n")
    
  })
  
  if (!is.null(outDir)) invisible(dev.off())
  
  if (returnResults){
    return(ret)
  }
  
  
  
}
