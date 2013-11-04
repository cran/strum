\name{strumVirtualModel-class}
\docType{class}

\alias{strumVirtualModel-class}
\alias{allRandomEffects<-}
\alias{allRandomEffects}
\alias{ascertainment<-}
\alias{ascertainment}
\alias{covariatePhenotypes}
\alias{formulas<-}
\alias{formulas}
\alias{paramNames<-}
\alias{paramNames}
\alias{varList<-}
\alias{varList}
\alias{allRandomEffects<-,strumVirtualModel-method}
\alias{allRandomEffects,strumVirtualModel-method}
\alias{ascertainment<-,strumVirtualModel-method}
\alias{ascertainment,strumVirtualModel-method}
\alias{covariatePhenotypes,strumVirtualModel-method}
\alias{formulas<-,strumVirtualModel-method}
\alias{formulas,strumVirtualModel-method}
\alias{paramNames<-,strumVirtualModel-method}
\alias{paramNames,strumVirtualModel-method}
\alias{varList<-,strumVirtualModel-method}
\alias{varList,strumVirtualModel-method}
\alias{show,strumVirtualModel-method}

\title{Class \code{"strumVirtualModel"}}
\description{
The \code{strumVirtualModel} class is a virtual class to represent a STRUctural 
Modeling of latent Variables for family data. It contains a description of the model as 
specified by the user through the model formulas and an set of functions to construct 
the internal matrix representation of the model.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
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
\section{Methods}{
  \describe{
    \item{formulas}{\code{signature(object = "strumVirtualModel")}: Accessor function, returns the formulas of the model.}
    \item{formulas<-}{\code{signature(object = "strumVirtualModel", value)}: Modifier function, sets the formulas of the model with a given value.}
    \item{covariatePhenotypes}{\code{signature(object = "strumVirtualModel")}: Accessor function, returns the names of covariate phenotypes of the model.}
    \item{varList}{\code{signature(object = "strumVirtualModel")}: Accessor function, returns the variable list of the model as a data.frame.}
    \item{varList<-}{\code{signature(object = "strumVirtualModel", value)}: Modifier function, sets the variable list of the model with a given value.}
    \item{allRandomEffects}{\code{signature(object = "strumVirtualModel")}: Accessor function, returns the names of all random effects of the model.}
    \item{allRandomEffects<-}{\code{signature(object = "strumVirtualModel", value)}: Modifier function, sets the allRandomEffects of the model with a given value.}
    \item{paramNames}{\code{signature(object = "strumVirtualModel")}: Accessor function, returns the names of all parameters of the model.}
    \item{paramNames<-}{\code{signature(object = "strumVirtualModel", value)}: Modifier function, sets the paramNames of the model with a given value.}
    \item{ascertainment}{\code{signature(object = "strumVirtualModel")}: Accessor function, returns the ascertainment definition of the model.}
    \item{ascertainment<-}{\code{signature(object = "strumVirtualModel", value)}: Modifier function, sets the ascertainment of the model with a given value.}
    \item{show}{\code{signature(object = "strumVirtualModel")}: Print a summary of the model.}
  }
}
\seealso{
\code{\linkS4class{strumModel}},  \code{\linkS4class{strumSimModel}}
}
\examples{
showClass("strumVirtualModel")
}
\keyword{classes}
