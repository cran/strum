\name{strumSimModel-class}
\docType{class}

\alias{strumSimModel-class}
\alias{markerInfo<-}
\alias{markerInfo}
\alias{traitMissingRate<-}
\alias{traitMissingRate}
\alias{markerInfo<-,strumSimModel-method}
\alias{markerInfo,strumSimModel-method}
\alias{traitMissingRate<-,strumSimModel-method}
\alias{traitMissingRate,strumSimModel-method}
\alias{show,strumSimModel-method}
\alias{plot,strumSimModel-method}

\title{Class \code{"strumSimModel"}}
\description{
\code{strumSimModel} is an S4 class that contains the information of a trait model 
for simulation.
}
\section{Objects from the Class}{
Objects should not be created by calls of the form \code{new("strumSimModel", ...)} 
but by the calls to the \code{\link{createSimModel}} function. 
}
\section{Slots}{
  \describe{
    \item{\code{markerInfo}:}{Object of class \code{"ANY"} containing either a strumMarker object or NULL.}
    \item{\code{traitMissingRate}:}{Object of class \code{"numeric"} containing the missing rate(s) of the simulated traits in the model.}
    \item{\code{varList}:}{Object of class \code{"data.frame"} containing a list of variables and their properties.}
    \item{\code{formulas}:}{Object of class \code{"character"} containing the information about the formulas and dependencies.}
    \item{\code{allRandomEffects}:}{Object of class \code{"character"} listing all variance components in the model.}
    \item{\code{paramNames}:}{Object of class \code{"character"} describing the names of the model parameters.}
    \item{\code{ascertainment}:}{Object of class \code{"ANY"} stating the ascertainment criteria (function) of the data.}
    \item{\code{E}:}{Object of class \code{"list"} containing the variance matrices of measurement errors of the model.}
    \item{\code{Z}:}{Object of class \code{"list"} containing the variance matrices of random effects for the structural equation of the model.}
    \item{\code{L}:}{Object of class \code{"function"} to generate a matrix to relate the unobserved latent factors to the observed traits in the model.}
    \item{\code{B}:}{Object of class \code{"function"} to generate a matrix denoting the causal relationship among the latent factors in the model.}
    \item{\code{Gs}:}{Object of class \code{"function"} to generate a matrix to relate the observed covariates to the unobserved latent factors in the model.}
    \item{\code{Gm}:}{Object of class \code{"function"} to generate a matrix to relate the observed covariates to the observed traits in the model.}
    \item{\code{thToThB}:}{Object of class \code{"numeric"} for future use.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{strumVirtualModel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{markerInfo}{\code{signature(object = "strumSimModel")}: Accessor function, returns the markerInfo of the model.}
    \item{traitMissingRate}{\code{signature(object = "strumSimModel")}: Accessor function, returns the traitMissingRate of the model.}
    \item{show}{\code{signature(object = "strumSimModel")}: Print a summary of the model.}
    \item{plot}{\preformatted{signature(object = "strumSimModel", layoutType="dot",
    name="strumSimModel", toFile=TRUE, fileType="dot")}: Plot the model.}
  }
}
\seealso{
\code{\linkS4class{strumModel}},  \code{\linkS4class{strumVirtualModel}},
\code{\link{createSimModel}}
}
\examples{
showClass("strumSimModel")
}
\keyword{classes}
