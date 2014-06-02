\name{createStrumModel}
\alias{createStrumModel}
\title{
 Create strumModel Object
}
\description{
A function to create an object of class strumModel.
}
\usage{
createStrumModel(formulas,
                 ascertainment = NULL,
                 defaultError = '<p,e>',
                 assumeExogCovariate = TRUE,
                 fixLoadingToOne = TRUE)
}
\arguments{
  \item{formulas}{Vector of expressions that defines the relationship among the 
  variables in the model (see Details).}
  \item{ascertainment}{Charactor stating the name of the column in data file that 
  contains the indicator variable (1, 0) designating the probands of the pedigrees.}
  \item{defaultError}{Vector of chractor listing the default variance components 
  in the model (see Details).}
  \item{assumeExogCovariate}{Logical. If TRUE, then it is assumed that any observed, 
  quantitative exogenous traits are covariates.  Otherwise, they are treated as  
  latent and dependent variables.}
  \item{fixLoadingToOne}{Logical. If TRUE, then it sets the first indicator variable 
  for each factor to have a coefficient of 1.  Otherwise, all coefficients are estimated.}  
}
\details{
This function is used to create a trait model which will be an input parameter 
for \code{strum} function.  

The formulas argument is given as a charactor string of several expressions that 
defines the relationship among the variables.  Blank lines and comments (line 
start with #) can be included between formulas.

Three different equations are allowed in the formulas syntax:

\enumerate{
  \item Measurement equations

  The "=~" operator may be read as "measured by" and it specifies the measurement 
  model equations in the model.  The left hand side of the "=~" must contain one 
  unobserved or latent factors.  The right hand side of the "=~" are observed 
  variables, observed covariates, and/or any measurement errors enclosed with "<>".
  
  \item Structural equations

  The "~" operator specifies the structural equations in the model.  The left hand 
  side of the "~" must contain one unobserved or latent factor variable name.  
  The right hand side of the "~" are unobserved or latent factors, observed 
  covariates, and/or any random effects enclosed with "<>".

  \item Adding a covariance

  By default, the endogenous error terms are uncorrelated, but the exogenous 
  variables are correlated.  To add a correlation between endogenous terms, use 
  "cov" operator.  For example, \preformatted{cov(x,y) = NA}.

  \item Removing a constraint on a coefficient

  By default, the value of \code{fixLoadingToOne} equals TRUE, setting 
  the first indicator variable for each factor to have a coefficient of 1.  To 
  disable this option selectively for a coefficient to be estimated, use "coef" 
  operator.  For example, \preformatted{coef(y1,z1) = NA}.

  \item Constraints

  The "=" operator specifies the constraints in the model, i.e., fixing a model 
  parameter - fixed variance, fixed covariance, or fixed coefficient.  
  The left hand side of the "=" must contain a reserved word for constraint - 
  "var", "cov", or "coef" - with one or two variable  names depending on the word.  
  The right hand side of the "=" is a fixed value.  Here are some examples:

  \itemize{

    \item Fixing the variance of a latent variable: 
      \preformatted{var(stress) = .1}

    \item Fixing the covariance between two latent variables:
      \preformatted{cov(z1,z2) = 4}

    \item Fixing the coefficient between two variables:
      \preformatted{coef(y1,z1) = 2}
  }
}

If a random effect such as \code{p}, \code{e}, \code{a} or \code{c} 
is not included in the model formulas, then, by default, the value of 
\code{defaultError} is included in the model.  For quantitative traits, the 
\code{e} term should always be present.  Therefore, the program automatically 
include \code{e} term even if no error terms are specified either in 
\code{formulas} nor in \code{defaultError}.

Note again that "a" (additive), "p" (polygenic), "c" (common environmental) and 
"e"(independent environmental) are reserved variable names to specify a type of 
variance component, so that may not be used as input arguments.
}
\value{
Returns an object of type \code{strumModel}.
}
\seealso{
\code{\linkS4class{strumModel}}
}
\examples{
# Model formulas.
#----------------
strumForm = 'bp     =~ SBP + DBP
             anger  =~ A1 + A2 + A3
             stress =~ S1 + S2 + S3
             bp ~ anger + stress + <p,e>
             stress ~ anger + rs6040343   
             var(stress)=.1  
            '

# Create a strumModel.
#---------------------
myStrumModel = createStrumModel(
                 formulas = strumForm,
                 ascertainment = "disease"
               )
}

