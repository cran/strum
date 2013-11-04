\name{strum}
\alias{strum}
\title{
STRUctural Modeling of Latent Variables for Family Data
}
\description{
Main function to run strum analysis and return the result as a (list of) strumFittedModel.
}
\usage{
strum(myStrumModel, myStrumData, ibdMarkers=NULL)
}
\arguments{
  \item{myStrumModel}{Object of class \code{strumModel}.}
  \item{myStrumData}{Object of class \code{strumData}.}
  \item{ibdMarkers}{Character vector contaning the name of marker(s) for additive genetic variance component \code{a} (see Details).}
}
\details{
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

