#==============================================================================
# File: strumVirtualModel.R
#
# Author: Nathan Morris
#
# Notes: strumVirtualModel class definition & methods
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Definition of strumVirtualModel class
#------------------------------------------------------------------------------

setClass("strumVirtualModel",
         representation(
           varList          = "data.frame",
           formulas         = "character",
           allRandomEffects = "character",
           paramNames       = "character",
           ascertainment    = "ANY",
           E                = "list",
           Z                = "list",
           L                = "function",
           B                = "function",
           Gs               = "function",
           Gm               = "function",
           thToThB          = "numeric"),
         representation("VIRTUAL"),
         prototype = list(
           varList          = data.frame(), #names(varList)=c("name","type","obs","covariate","exogen")
           formulas         = character(),
           allRandomEffects = character(),
           paramNames       = character(),
           ascertainment    = NULL,
           E                = list(),
           Z                = list(),
           L                = function(onefamilydf) return(TRUE),
           B                = function(onefamilydf) return(TRUE),
           Gs               = function(onefamilydf) return(TRUE),
           Gm               = function(onefamilydf) return(TRUE),
           thToThB          = numeric())
) 

#------------------------------------------------------------------------------
# 'covariatePhenotypes' accessor functions:
#------------------------------------------------------------------------------
setGeneric('covariatePhenotypes', function(object) standardGeneric('covariatePhenotypes'))
setMethod('covariatePhenotypes', 'strumVirtualModel',
          function(object)
          {
            return(as.character(object@varList$name[object@varList$covariate & object@varList$distribution!="mendelian"]))
          }
)

#------------------------------------------------------------------------------
# 'varList' accessor functions:
#------------------------------------------------------------------------------
setGeneric('varList', function(object) standardGeneric('varList'))
setMethod('varList', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@varList)
          }
)

setGeneric('varList<-', function(object,value) standardGeneric('varList<-'))
setMethod('varList<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@varList <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'formulas' accessor functions:
#------------------------------------------------------------------------------
setGeneric('formulas', function(object) standardGeneric('formulas'))
setMethod('formulas', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@formulas)
          }
)

setGeneric('formulas<-', function(object,value) standardGeneric('formulas<-'))
setMethod('formulas<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@formulas <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'allRandomEffects' accessor functions:
#------------------------------------------------------------------------------
setGeneric('allRandomEffects', function(object) standardGeneric('allRandomEffects'))
setMethod('allRandomEffects', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@allRandomEffects)
          }
)

setGeneric('allRandomEffects<-', function(object,value) standardGeneric('allRandomEffects<-'))
setMethod('allRandomEffects<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@allRandomEffects <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'paramNames' accessor functions:
#------------------------------------------------------------------------------
setGeneric('paramNames', function(object) standardGeneric('paramNames'))
setMethod('paramNames', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@paramNames)
          }
)

setGeneric('paramNames<-', function(object,value) standardGeneric('paramNames<-'))
setMethod('paramNames<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@paramNames <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# 'ascertainment' accessor functions:
#------------------------------------------------------------------------------
setGeneric('ascertainment', function(object) standardGeneric('ascertainment'))
setMethod('ascertainment', signature(object = 'strumVirtualModel'),
          function(object)
          {
            return(object@ascertainment)
          }
)

setGeneric('ascertainment<-', function(object,value) standardGeneric('ascertainment<-'))
setMethod('ascertainment<-', signature(object = 'strumVirtualModel'),
          function(object, value)
          {
            object@ascertainment <- value
            return(object)
          }
)

#------------------------------------------------------------------------------
# show generic functions
#------------------------------------------------------------------------------
setMethod("show", signature(object = "strumVirtualModel"),
          function(object) 
          {
            .showModel(object, "strumVirtualModel")
          }
)

#------------------------------------------------------------------------------
# Common function to show Model class 
#------------------------------------------------------------------------------
.showModel = function(object, name) 
{
  cat("\nBasic properties of the model:\n")
    .printInfoLine("Model Class", name, 40)

  if( is.null(object@ascertainment) )
    .printInfoLine("Ascertainment", "FALSE", 40)
  else
    .printInfoLine("Ascertainment", "TRUE", 40)

  cat("\nList of all variables:\n")
  myvarList = object@varList
  myrowNames = myvarList$name
  myvarList = myvarList[,names(myvarList)!="endogenous" & names(myvarList)!="name" ]
  mynames = names(myvarList)
  mynames[mynames=="variableLevels"] = "Level"
  mynames[mynames=="exogenous"]      = "Exogen."
  mynames[mynames=="distribution"]   = "Distr."
  mynames[mynames=="parameter"]      = "Param."
  mynames = paste(toupper(substring(mynames, 1,1)), substring(mynames, 2), sep="")
  names(myvarList) = mynames
  print(myvarList, row.names = myrowNames)

  cat("\nModel formulas:\n")
  cat(object@formulas)
}
