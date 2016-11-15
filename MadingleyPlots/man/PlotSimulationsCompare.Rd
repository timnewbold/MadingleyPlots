\name{PlotSimulationsCompare}
\alias{PlotSimulationsCompare}
\title{Plot Impacts from Madingley Simulations}
\description{
Plots properties of ecosystems under impact, compared with a specified baseline simulation
}
\usage{
  PlotSimulationsCompare <- function(resultsDirBaseline,resultsDir,plotName,outDir,
                               labelBaseline,labels,xVals,vars,cols,ylab,xlab,
                               confidenceInterval,endTimeStepBaseline,
                               numTimeStepsBaseline,endTimeStep,numTimeSteps,
                               whichCells,returnResults,pairedSimulations)
}
\arguments{
  \item{resultsDirBaseline}{
    The directory containing the baseline simulation
  }
  \item{resultsDir}{
    The directory containing the impact simulations (can be the same as the baseline directory)
  }
  \item{plotName}{
    The name to be given to the output file
  }
  \item{outDir}{
    The directory to which to save the output file
  }
  \item{labelBaseline}{
    The simulation label identifying the baseline simulations. The default is null, in which case the function will attempt to identify the label based on the filenames
  }
  \item{labels}{
    A vector of simulation labels identifying the impact simulations. The default is null, in which case the funciton will attempt to identify a single label based on the filenames
  }
  \item{xVals}{
    The impact values for each impact level (used as x-axis values). The default is null, in which case values will be from 1 to the number of impact levels
  }
  \item{vars}{
    The names of the ecosystem varaibles to plot. Defaults are the biomasses for each trophic level.
  }
  \item{cols}{
    The colours to use for the error bars. Defaults are the default colours for the default variables
  }
  \item{ylab}{
    The y-axis label to use. The default is to be the same as plotName
  }
  \item{xlab}{
    The x-axis label to use. Default is to have no x-axis label
  }
  \item{confidenceInterval}{
    The confidence interval to use for the error bars. Default is 95%
  }
  \item{endTimeStepBaseline}{
    The final time step from the baseline simulations to consider in the calculations. Default is the final time step in the simulation
  }
  \item{numTimeStepsBaseline}{
    The number of time steps from the baseline simulations to consider in the calculations. The default is 12 (1 year assuming a monthly time step)
  }
  \item{endTimeStep}{
    The final time step from the impact simulations to consider in the calculations. Default is the final time step in the simulation
  }
  \item{numTimeSteps}{
    The number of time steps from the impact simulations to consider in the calculations. The default is 12 (1 year assuming a monthly time step)
  }
  \item{whichCells}{
    A vector indicating which cells within the simulations to plot. Default is to plot all cells
  }
  \item{returnResults}{
    Whether to return a data frame containing the results
  }
  \item{pairedSimulations}{
    Whether the baseline and impact simulations are paired with one another (for example if the same initialisation simulations were used for both). Default (and the only option supported at present) is for paired simulations
  }
}
\details{
  \author{
    Tim Newbold <t.newbold@ucl.ac.uk>
  }
}