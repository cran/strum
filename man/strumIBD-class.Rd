\name{strumIBD-class}
\docType{class}

\alias{strumIBD-class}
\alias{ibdMarker<-}
\alias{ibdMarker}
\alias{ibdMatrix<-}
\alias{ibdMatrix}
\alias{ibdMarker<-,strumIBD-method}
\alias{ibdMarker,strumIBD-method}
\alias{ibdMatrix<-,strumIBD-method}
\alias{ibdMatrix,strumIBD-method}
\alias{show,strumIBD-method}

\title{Class \code{"strumIBD"}}
\description{
\code{strumIBD} is an S4 class that represents IBD information of the input data for 
a STRUctural Modeling of latent Variables for family data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("strumIBD", ...)} 
or by the \code{importIBD} method in \code{\linkS4class{strumData}} class. 
}
\section{Slots}{
  \describe{
    \item{\code{ibdMarker}:}{Object of class \code{"data.frame"} containing the information about the markers.  Typically contains two columns, "Marker" and "Position".}
    \item{\code{ibdMatrix}:}{Object of class \code{"list"} containing the lists of identical by descent (ibd) matrices, for each family and marker.}
  }
}
\section{Methods}{
  \describe{
    \item{ibdMarker}{\code{signature(object = "strumIBD")}: Accessor function, returns the marker information.}
    \item{ibdMatrix}{\code{signature(object = "strumIBD")}: Accessor function, returns the list of ibd matrices of the input data.}
    \item{show}{\code{signature(object = "strumIBD")}: Print a summary of the class.}
  }
}
\seealso{
\code{\linkS4class{strumSimModel}},  \code{\linkS4class{strumData}}
}
\examples{
showClass("strumIBD")
}
\keyword{classes}
