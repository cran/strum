#==============================================================================
# File: strum.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Main driver of strum analysis.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# The main analysis given StrumModel and input data
#------------------------------------------------------------------------------
strum = function(myStrumModel, myStrumData, ibdMarkers=NULL)
{
  # 1. Check input parameters
  #---------------------------
  if( class(myStrumModel) != "strumModel" )
    stop("Wrong class input! Need an object of 'strumModel' class.")

  if( class(myStrumData) != "strumData" )
    stop("Wrong class input! Need an object of 'strumData' class.")

  isPedigree = (dataType(myStrumData) == "Pedigree")

  # 2. Check ascertainment
  #------------------------
  nonProbands = list()

  if( isPedigree )
  {
    aName = ascertainment(myStrumModel)

    if( class(aName) == "character" )
    {
      if( !(aName %in% names(dataVals(myStrumData))) )
      {
        stop(paste("Model contains ascertainment but no column of that name exist in data. ",
                   "Either set ascertainment in the model to NULL",
                   "or use a data set with the required proband column.",
                   sep = " "))
      } else
      {
        aValues = dataVals(myStrumData)[,aName, drop=FALSE]
        aValues = split(aValues, dataVals(myStrumData)$family)

        nonProbands = lapply(aValues, function(pk) return(which(pk==0)))

        pCount = unlist(lapply(aValues, sum))
        if( any(pCount > 1) )
          warning("Pedigree(s) with more than 1 proband exist!")
        if( any(pCount== 0) )
          warning("Pedigree(s) with no proband exist!")
      }

    } else if( !is.null(aName) )
    {
      warning("Model contains a wrongly specified ascertainment!  Ignoring...")
    }
  }

  # 3. Check random effect
  #------------------------
  allRE = allRandomEffects(myStrumModel)

  if( any(c("a","p","c") %in% allRE) & !isPedigree )
    stop("Can't fit a genetic, pologenic or common environment effect without pedigree structures.")

  cat("\n Start STRUM analysis ...\n")

  # 3. Prepare analysis data
  #--------------------------
  y     = .getAnalysisY( myStrumModel, myStrumData, nonProbands)
  x     = .getAnalysisX( myStrumModel, myStrumData, nonProbands)
  vcPEC = .getAnalysisVC(myStrumModel, myStrumData, allRE)

  myFittedModel = list()

  if( any(allRE == "a") )
  {
    iMarkers = NULL
    
    if( is.null(ibdMarkers) )
    {
      iMarkers = names(ibdMatrix(ibd(myStrumData)))
    } else
    {
      for( m in 1:length(ibdMarkers) )
      {
        mName = ibdMarkers[m]

        if( !(mName %in% names(ibdMatrix(ibd(myStrumData)))) )
          warning(paste("IBD marker \"", mName, "\" doesn't exist in the IBD file!  Ignoring...",
                        sep=""))
        else
          iMarkers = c(iMarkers, mName)
      }

      if( is.null(iMarkers) )
      {
        stop(paste("No valid IBD marker available in 'ibdMarkers'!  ",
                   "Will use all markers in 'myStrumData' object.",
                   sep=""))

        iMarkers = names(ibdMatrix(ibd(myStrumData)))
      }
    }

    if( is.null(iMarkers) )
    {
      stop(paste("No IBD marker available in 'myStrumData' object!  Need an object of ",
                 "strumData class constructed with the IBD information of the marker(s) ",
                 "to be included in the model as additive genetic component(s).",
                 sep=""))
    }

    for( m in 1:length(iMarkers) )
    {
      mName = iMarkers[m]

      cat(paste("\n  ", mName, ":\n", sep=""))
      vca = ibdMatrix(ibd(myStrumData))[[mName]]

      vcAll = mapply(function(peci, vcai)
                     {
                       peci[["a"]] = vcai
                       return(peci)
                     },
                     vcPEC, vca,
                     SIMPLIFY = FALSE)

      vcPro = list()
      if( length(nonProbands) > 0 )
        vcPro = mapply(.filterMissingVC, vcPEC, nonProbands, nonProbands, SIMPLIFY=FALSE) 

      vc = list(vcAll = vcAll, vcPro = vcPro)
        
      mFittedModel = .fitModel(myStrumModel, y, x, vc)

      myFittedModel[[mName]] = mFittedModel
    }
  } else
  {
    vcPro = list()
    if( length(nonProbands) > 0 )
      vcPro = mapply(.filterMissingVC, vcPEC, nonProbands, nonProbands, SIMPLIFY=FALSE) 

    vc = list(vcAll = vcPEC, vcPro = vcPro)

    myFittedModel = .fitModel(myStrumModel, y, x, vc)
  }

  cat("\nAnalysis completed!\n\n")

  return(myFittedModel)
}

#------------------------------------------------------------------------------
# Model analysis data Y from StrumData given StrumModel
# - Extract Y
#------------------------------------------------------------------------------
.getAnalysisY = function(myStrumModel, myStrumData, nonProbands)
{
  yNames = myStrumModel@varList$name[myStrumModel@varList$inY==TRUE]

  # Error check if names not there
  if( !length(yNames) )
    stop("No y names exist in StrumModel")

  yAll = dataVals(myStrumData)[,yNames, drop=FALSE]
  yAll = split(yAll, dataVals(myStrumData)$family)
  yAll = lapply(yAll, data.matrix)
  
  yPro = list()
  if( length(nonProbands) > 0 )
    yPro = mapply(.filterMissing, yAll, nonProbands, TRUE, SIMPLIFY=FALSE) 

  return(list(yAll = yAll, yPro = yPro))
}

#------------------------------------------------------------------------------
# Model analysis data X from StrumData given StrumModel
# - Extract X from covariates
#------------------------------------------------------------------------------
.getAnalysisX = function(myStrumModel, myStrumData, nonProbands)
{
  xNames = myStrumModel@varList$name[myStrumModel@varList$covariate==TRUE]

  xAll = list()

  # If names not there, then intercept only
  if( !length(xNames) )
  {
    xAll = rep(1, nrow(dataVals(myStrumData)))
    xAll = split(xAll, dataVals(myStrumData)$family)
    xAll = lapply(xAll, data.matrix)
  } else
  {
    xAll = dataVals(myStrumData)[,xNames, drop=FALSE]
    xAll = lapply(split(xAll, dataVals(myStrumData)$family),
                  function(xi)
                  {
                    xi = data.matrix(xi)
                    model.matrix(~xi)
                  })
  }
  
  xPro = list()
  if( length(nonProbands) > 0 )
    xPro = mapply(.filterMissing, xAll, nonProbands, TRUE, SIMPLIFY=FALSE) 

  return(list(xAll = xAll, xPro = xPro))
}

#------------------------------------------------------------------------------
# Model analysis data VC from StrumData given StrumModel
# - Use allRandomEffects to extract VC values
#------------------------------------------------------------------------------
.getAnalysisVC = function(myStrumModel, myStrumData, allRE)
{
  vcPEC = list()

  if( any(allRE == "p") )
    vcPEC = lapply(phi(myStrumData), function(p) return(list(p=p)))

  if( any(allRE == "e") )
  {
    vce = lapply(split(dataVals(myStrumData), dataVals(myStrumData)$family),
                 function(p) return(list(e=diag(nrow(p)))))

    if( length(vcPEC) )
      vcPEC = mapply(function(vci, vcei) return(c(vci, vcei)),
                     vcPEC, vce,
                     SIMPLIFY = FALSE)
    else
      vcPEC = vce
  }

  if( any(allRE == "c") )
  {
    vcc = lapply(phi(myStrumData), 
                 function(p) return(list(c=matrix(rep(1,nrow(p)*nrow(p)),nrow=nrow(p)))))

    if( length(vcPEC) )
      vcPEC = mapply(function(vci, vcci) return(c(vci, vcci)),
                     vcPEC, vcc,
                     SIMPLIFY = FALSE)

    else
      vcPEC = vcc
  }

  return(vcPEC)
}
