\name{strum}
\alias{strum}
\title{
STRUctural Modeling of Latent Variables for Family Data
}
\description{
Main function to run strum analysis and return the result as a (list of) strumFittedModel.
}
\usage{
strum(myStrumModel,
      myStrumData,
      step1OptimControl = list(maxit=5500, fnscale=-10),
      startValueControl = list(initPopulation=NULL,
                               nChildren=NULL,
                               nGenerations=NULL,
                               selection1=NULL,
                               selection2=NULL),
      step2OptimControl = list(maxit=5000, reltol=.Machine$double.eps),
      ibdMarkers=NULL)
}
\arguments{
  \item{myStrumModel}{Object of class \code{strumModel}.}
  \item{myStrumData}{Object of class \code{strumData}.}
  \item{step1OptimControl}{List of control parameters for \code{optim} function in step1.}
  \item{startValueControl}{List of control parameters for strating value generation function in step2.}
  \item{step2OptimControl}{List of control parameters for \code{optim} function in step2.}
  \item{ibdMarkers}{Character vector contaning the name of marker(s) for additive genetic variance component \code{a} (see Details).}
}
\details{
Three parameters of this function provide more flexibility to users to control the 
strum analysis in different stage.

The \code{step1OptimControl} argument is given as a list of several control parameters 
that is passed directly to the function \code{optim} in the step 1 of model fitting. 

The \code{step2OptimControl} argument is given as a list of several control parameters  
that is passed directly to the function \code{optim} in the step 2 of model fitting. 

The \code{startValueControl} argument is given as a list of 5 parameters for the function 
to generate the starting values (using genetic algorithm) of model fitting in step 2.  Five different parameters in the list are:

\itemize{
  \item \code{initPopulation}

  The size of initial population. The default value is the number of model parameters * 120.

  \item \code{nChildren}

  The number of children. The default value is the number of model parameters * 2.

  \item \code{nGenerations}

  The number of generation. The default value is 20.

  \item \code{selection1}

  The size of the initial selection. The default value is 15.

  \item \code{selection2}

  The size of the final selection. The default value is 15.
}

Note that if the analysis model includes an additive genetic variance component \code{a} 
for the random effects, the value of \code{iMarkers} is used to determine the name of 
marker(s) for \code{a} from the imported ibd file.  If no values are specified for 
\code{ibdMarkers}, then, by default, all imported ibd markers that exist in the  
\code{myStrumData} object are used, analyzing one by one.  Therefore, the analysis result 
includes a list of strumFittedModel.
}
\value{
Returns a (list of) strumFittedModel object.
}
\seealso{
\code{\link{createStrumModel}}, \code{\link{createStrumData}}
}
\examples{
\dontrun{
# Run strum analysis.
# - myStrumModel is a strumModel object.
# - myStrumData is a strumData object.
#---------------------------------------
fitResult = strum(myStrumModel, myStrumData, ibdMarkers=c("marker1","marker2"))
}
}

