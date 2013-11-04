\name{strumData-class}
\docType{class}

\alias{strumData-class}
\alias{dataType}
\alias{dataType<-}
\alias{dataVals}
\alias{dataVals<-}
\alias{ibd}
\alias{phi}
\alias{phi<-}
\alias{importIBD}
\alias{dataType,strumData-method}
\alias{dataType<-,strumData-method}
\alias{dataVals,strumData-method}
\alias{dataVals<-,strumData-method}
\alias{ibd,strumData-method}
\alias{phi,strumData-method}
\alias{phi<-,strumData-method}
\alias{$,strumData-method}
\alias{$<-,strumData-method}
\alias{[,strumData,ANY,ANY,ANY-method}
\alias{[<-,strumData,ANY,ANY,ANY-method}
\alias{[[,strumData-method}
\alias{[[<-,strumData-method}
\alias{importIBD,strumData-method}
\alias{show,strumData-method}

\title{Class \code{"strumData"}}
\description{
\code{strumData} is an S4 class that contains the input data for a STRUctural Modeling 
of latent variables for family data.
}
\section{Objects from the Class}{
Objects should not be created by calls of the form \code{new("strumData", ...)} 
but by the calls to the \code{\link{createStrumData}} function. 
}
\section{Slots}{
  \describe{
    \item{\code{dataType}:}{Character stating the type of input data, must be 
    \code{"Pedigrees"} or \code{"RawData"}}
    \item{\code{dataVals}:}{Object of class \code{data.frame} containing input data.}
    \item{\code{phi}:}{Object of class \code{"list"} containing the kinship matrices for each pedigree.  This is empty if the dataType is not "Pedigree".}
    \item{\code{ibd}:}{Object of class \code{"strumIBD"} containing the IBD information for input pedigrees.  This is empty if the dataType is not "Pedigree".}
  }
}
\section{Methods}{
  \describe{
    \item{importIBD}{\code{signature(object = "strumData", fileName, fileType="SAGE")}: Import the ibd information of the input data.}
    \item{show}{\code{signature(object = "strumData")}: Print a summary of the class.}
  }
}
\seealso{
\code{\linkS4class{strumIBD}}, \code{\link{createStrumData}}, \code{\link{simulateStrumData}}
}
\examples{
showClass("strumData")
}
\keyword{classes}
