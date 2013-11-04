\name{createSimModel}
\alias{createSimModel}
\title{
Create strumSimModel Object
}
\description{
A function to create an object of class strumSimModel which contains the information 
used for simulation.
}
\usage{
createSimModel(formulas,
               tMissingRate = c(0),
               markerInfo = NULL,
               ascertainment = NULL,
               defaultError = '<p,e>',
               assumeExogCovariate = TRUE)
}
\arguments{
  \item{formulas}{Vector of expressions that defines the relationship among the 
  variables in the model with the fixed coefficient (see Details).}
  \item{markerInfo}{Object containing a strumMarker or NULL.}
  \item{tMissingRate}{Numeric vector to hold the missing rate(s) of the simulated 
  trait(s) in the model.  The length of vector is equal to the number of observed 
  traits or 1 if one missing rate is applied to all traits.}
  \item{ascertainment}{Function stating the ascertainment criteria of the data (see Details).}  
  \item{defaultError}{Vector of chractor listing the default variance components 
  in the model (see Details).}
  \item{assumeExogCovariate}{Logical. If TRUE, then it is assumed that any observed, 
  quantitative exogenous traits are covariates.  Otherwise, they are treated as  
  latent and dependent variables.}
}
\details{
This function is used to create a \code{strumSimModel} object which will be an 
input parameter for \code{simulateStrumData} function. 

The \code{formulas} argument is given as a charactor string of several expressions that 
defines the relationship among the variables.  Blank lines and comments (line start 
with #) can be included between formulas.

Note that the coefficients of the right hand side variables of the measurement 
equations are fixed with the value given in the formulas.

Three different equations are allowed in the formulas syntax:

\enumerate{
  \item Measurement equations

  The "=~" operator may be read as "measured by" and it specifies the measurement 
  equations in the model.  The left hand side of the "=~" must contain one 
  unobserved or latent factor variable name.  The right hand side of the "=~" are 
  observed variables, observed covariates, and/or any measurement errors enclosed 
  with "<>".

  \item Structural equations

  The "~" operator specifies the structural equations in the model.  The left hand 
  side of the "~" must contain one unobserved or latent factor variable name.  
  The right hand side of the "~" are unobserved or latent factors, observed 
  covariates, and/or any random effects enclosed with "<>".

  \item Constraints

  The "=" operator specifies the constraints in the model, i.e., fixing a model 
  parameter - fixed variance, fixed covariance.  
  The left hand side of the "~" must contain a reserved word for constraint - 
  "var" or "cov" - with one or two variable  names depending on the word.  
  The right hand side of the "=" is a fixed value.  Here are some examples:

  \itemize{

    \item Fixing the variance of a latent variable: 
      \preformatted{var(stress) = .1}

    \item Fixing the covariance between two latent variables:
      \preformatted{cov(z1,z2) = 4}
  }
}

If a random effect such as \code{p}, \code{e}, \code{a} or \code{c} 
is not included in the model formulas, then, by default, the value of 
\code{defaultError} is included in the model.  For quantitative traits, the 
\code{e} term should always be present.  Therefore, the program automatically 
include \code{e} term even if no error terms are specified either in 
\code{formulas} nor in \code{defaultError}.

Note that "a" (additive), "p" (polygenic), "c" (common environmental) and 
"e"(independent environmental) are reserved variable names to specify a type of 
variance component, so that may not be used as input arguments.

The return value of \code{ascertainment} function is composed of two components for 
each pedigree, allowing either one of components or both components to be returned as a list.  
The first component is a TRUE/FALSE value indicating the ascertainment status of the pedigree.  
The second component is a vector of TRUE/FALSE stating the proband status of each 
member of the pedigree.  Here are some examples:

\enumerate{
  \item Ascertainment status only

  \preformatted{aFunction = function(ped) return(TRUE)
aFunction = function(ped) return(any(ped$SBP > 2))}

  The pedigree is ascertained when a "TRUE" is returned for that pedigree.

  \item Proband status only

  \preformatted{aFunction = function(ped) return(ped$SBP > 2) 
aFunction = function(ped) return(ped$SBP > 2 && ped$DBP < 1)
aFunction = function(ped) return(ped$disease == 1)}

  In this case, the ascertainment status of the pedigree is determined according to the 
  existance of a proband in each pedigree.  If no probands who meet the criteria exist in 
  a simulated pedigree, then that pedigree is not selected.

  \item Both

  \preformatted{aFunction = function(ped)
{
  asc = any(ped$disease==1)    
  pro = ped$SBP > 2
  return(list(aStatus=asc, pStatus=pro))
}}

  When a list with both information is returned, it is required that the ascertainment 
  status (a logical value) is the first component, and the proband status (a vector of 
  logicals of pedigree size) is the second component.  The pedigree is selected when the 
  returned first component is a "TRUE" for that pedigree.
}
}
\value{
Returns an object of type \code{strumSimModel}.
}
\seealso{
\code{\linkS4class{strumSimModel}}, \code{\linkS4class{strumMarker}}
}
\examples{
\dontrun{
# Create a strumSimModel.
# - simForm is a formulas string.
# - snpStrumMarker is an object of strumMarker class.
# - aFunction is an ascertainment function definition.
#------------------------------------------------------
mySimModel = createSimModel(
                 formulas = simForm,
                 markerInfo = snpStrumMarker,
                 ascertainment = aFunction
               )
}
}

