\name{simulateStrumData}
\alias{simulateStrumData}
\title{
 Simulate strumData Object
}
\description{
A function to create a strumData object containing the simulated data according 
to the simulation model.
}
\usage{
simulateStrumData(simModel, inData=NULL, N=NULL)
}
\arguments{
  \item{simModel}{Object of class \code{strumSimModel}.}
  \item{inData}{Object of class \code{strumData}, \code{data.frame}, or NULL.}
  \item{N}{A positive interger number to specify the particular size of simulated data or NULL.}
}
\details{
This function is used to simulate a data according to the specified simulation model.  
Note that either inData or N, or both must be specified.
}
\value{
Returns an object of class \code{strumData}.
}
\seealso{
\code{\linkS4class{strumSimModel}}, \code{\linkS4class{strumData}}
}
\examples{
\dontrun{
# Simulate a strumData object with different type of input data.
# - mySimModel is an object of strumSimModel class.
# - dF is a data.frame containing input data.
# - rawStrumData is an object of "RawData" type strumData class.
# - pedStrumData is an object of "Pedigree" type strumData class.
#-----------------------------------------------------------------
myStrumData = simulateStrumData(mySimModel, N=1000)
myStrumData = simulateStrumData(mySimModel, dF)
myStrumData = simulateStrumData(mySimModel, rawStrumData)
myStrumData = simulateStrumData(mySimModel, pedStrumData)
}
}

