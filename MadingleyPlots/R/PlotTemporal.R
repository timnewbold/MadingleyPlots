
PlotTemporal<-function(resultsDir,plotName,outDir=NULL,
                       label=NULL,
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
                       returnResults = FALSE,
                       whichCells = NULL){
  
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
  
  # Find the unique cells in these simulations
  cells.re<-regexpr("Cell[0-9]+",files)
  cells<-as.list(unique(substr(files,cells.re,cells.re+
                         attr(cells.re,"match.length")-1)))
  
  if (!is.null(whichCells)){
    cells <- cells[whichCells]
  }
  
  # Find the simulation numbers
  sims.re<-regexpr("_[0-9]+_",files)
  sims<-as.list(unique(substr(files,sims.re,sims.re+
                        attr(sims.re,"match.length")-1)))
  
  if(is.null(label)){
    label<-unique(substr(files,1,sims.re-1))
    stopifnot(length(label)==1)
    label<-label[1]
  } else {
    label <- paste("BasicOutputs_",label,sep="")
  }
  
  .Log(paste("Found results for ",length(cells)," cells\n",sep=""))
  .Log(paste("Found results for ",length(sims)," simulations\n",sep=""))
  
  .Log("Getting basic information about simulations\n")
  sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sims[1],cells[1],
                  ".nc",sep="")
  data<-open.sds(sds.path)
  times<-get.sds(data,"Time step")
  y<-times/12
  y<-ceiling(y)
  years<-unique(y)
  
  .Log("Initializing plot\n")
  dims<-.plotDims(length(cells))
  if(!is.null(outDir)){
    pdf(paste(outDir,plotName,".pdf",sep=""),
        width = dims$width,height = dims$height)
  }
    
  par(mfrow=.gridArrange(length(cells)))
  par(mar=c(3,3.5,0.5,0.5))
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
      sds.path<-paste("msds:nc?file=",resultsDir,"/",label,sim,cell,
                      ".nc",sep="")
      
      data<-open.sds(sds.path)
      
      # Populate the results matrices
      for (var in vars){
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
    plot(999999999999999,999999999999999,xlim=range(years),ylim=c(minVal,maxVal),log="y",
         xlab=xlab,yaxt="n",ylab=NA)
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
    
    
  })
  
  if (!is.null(outDir)) invisible(dev.off())
  
  if (returnResults){
    return(ret)
  }
  
  
  
}
