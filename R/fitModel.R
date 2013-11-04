#==============================================================================
# File: fitModel.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Main driver of fitting the model.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Fitting the Model.
#------------------------------------------------------------------------------
.fitModel = function(model, y, x, vc)
{  
  numTrait = ncol(y$yAll[[1]])
  numCov   = ncol(x$xAll[[1]])
  numVC    = length(vc$vcAll[[1]])
  
  numC = numCov * numTrait
  numV = numVC * (numTrait*(numTrait+1)/2)

  df0 = (numC+numV) - length(paramNames(model))

  # error - Not identifiable, more parameters in model than in saturated model!
  if( df0 < 0 )
  {
    myFittedModel = new("strumFittedModel",
                         myStrumModel  = model,
                         modelValidity = 1)

    return(myFittedModel)
  }

  # 1. Fit stage1 
  #---------------
  s1 = .fitModelStep1(y, x, vc)

  # 2. Fit stage2 
  #---------------
  s2 = .fitModelStep2(s1, model)

  # error - Parameters are not locally identifiable!
  if( df0 != ncol(s2$dfOrthog) )
  {
    myFittedModel = new("strumFittedModel",
                         myStrumModel  = model,
                         modelValidity = 2)

    return(myFittedModel)
  }

  # 3. Test the model
  #-------------------
  tOut = .testModel(model, s1, s2)

  fitVal = 0

  # Same number of parameters in model and saturated model.  No fit test!
  if( !is.data.frame(tOut$parEst) )
  {
    fitVal = 4
  } else if( !is.data.frame(tOut$fitStat) )
  {
    fitVal = 3
  }

  myFittedModel = new("strumFittedModel",
                       myStrumModel              = model,
                       modelValidity             = fitVal,
                       fittedParameters          = tOut$parEst,
                       fittedParametersCovMatrix = s2$thetaCov,
                       deltaParameters           = s1$delta,
                       deltaParametersCovMatrix  = s1$deltaCov,
                       parDiff                   = s2$diff,
                       parDiffCovMatrix          = s2$diffCov,
                       chiTestOut                = tOut$fitStat)

  return(myFittedModel)
}

#------------------------------------------------------------------------------
# Fitting the Model, stage 1
# - Form a limited information estimate for the saturated model.
#------------------------------------------------------------------------------
.fitModelStep1 = function(y, x, vc)
{
  CV = .estimateDeltaParameter(y, x, vc)
  De = .estimateDeltaDerivative(y, x, vc, CV$C, CV$V)
  JF = .estimateDeltaCovariance(De)

  delta = .Call("getDelta", t(CV$C), CV$V)

  .printInfoLine("Fitting model step 1", "Done", 50, 2)

  return(list(K=JF$K, delta=delta, deltaCov=JF$deltaCov, w=JF$w))
}

#------------------------------------------------------------------------------
# Fitting the Model, stage 2
# - Estimate parameters and asymptotic covariance matrix of the model.
#------------------------------------------------------------------------------
.fitModelStep2 = function(s1, model)
{
  startTheta = .generateStartValue(model, s1$delta, s1$w)
  optimTheta = .estimateThetaParameter(model, s1$delta, s1$w, startTheta)
  parCov     = .estimateThetaCovariance(model, s1, optimTheta)

  .printInfoLine("Fitting model step 2", "Done", 50, 2)

  return(list(theta    = optimTheta,
              thetaCov = parCov$thetaCov,
              diff     = parCov$diff,
              diffCov  = parCov$diffCov,
              dfOrthog = parCov$dfOrthog,
              U        = parCov$U))
}

#------------------------------------------------------------------------------
# Test the model fit.
#------------------------------------------------------------------------------
.testModel = function(model, s1, s2)
{
  # Calculate SE, CI, and p values of theta.
  #------------------------------------------------------------------------------
  theta    = s2$theta
  thetaCov = s2$thetaCov

  myeps = 1.0e-11 #.Machine$double.eps (2.220446e-16) or 1.0e-11

  thetaVar = diag(thetaCov)
  thetaVar[which(abs(thetaVar) < myeps)] = 0

  if( length(which(thetaVar < 0)) > 0 )
  {
    return(list(parEst  = NA,
                fitStat = NA))
  }
    
  standardError = sqrt(thetaVar)

  lowerCI = theta + qnorm(.025) * standardError
  upperCI = theta - qnorm(.025) * standardError

  z      = theta/standardError
  pValue = pnorm(-abs(z))

  dfPar = data.frame(estimate = theta,
                     stdError = standardError,
                     lowerCI  = lowerCI,
                     upperCI  = upperCI,
                     pValue   = pValue)

  if( length(s1$delta) == length(theta) )
  {
    return(list(parEst  = dfPar,
                fitStat = NA))
  }

  # Fit check - Chi-square statistics
  #-----------------------------------
  T0       = sum(s2$diff*s2$diff*s1$w) * s1$K
  df0      = length(s1$delta) - length(theta)
  dfOrthog = s2$dfOrthog
  gamma    = s1$deltaCov
  Ugamma   = s2$U %*% gamma
  trUgamma = sum(diag(Ugamma))

  # 1. Mean scaled chi-square statistics
  #--------------------------------------
  temp     = (t(dfOrthog) %*% gamma %*% dfOrthog)
  out      = qr(temp)
  rank     = out$rank

  kappa1 = trUgamma / rank
  T1     = T0 / kappa1

  # 2. Mean/variance scaled chi-square statistics
  #-----------------------------------------------
  dstar = trUgamma^2 / sum(diag(Ugamma^2))
  d2    = round(dstar)

  kappa2 = trUgamma / d2
  T2     = T0 / kappa2

  pVal0 = pchisq(T0, length(theta), lower.tail = FALSE)
  pVal1 = pchisq(T1, rank,          lower.tail = FALSE)
  pVal2 = pchisq(T2, d2,            lower.tail = FALSE)

  # 3. Corrected chi-square statistics
  #------------------------------------
  M = s2$diffCov

  eM = eigen(M, symmetric=TRUE)
  eM$values[which(abs(eM$values) < myeps)] = 0
  eM$values[which(eM$values < 0)] = 0 
  eM$vectors[which(abs(eM$vectors) < myeps)] = 0
  
  eMsqrt = eM$vectors %*% diag(sqrt(eM$values)) %*% t(eM$vectors)

  B = t(eMsqrt) %*% diag(s1$w) %*% eMsqrt

  eB = eigen(B, symmetric=TRUE)

  nEival = length(eB$values)
  nsims  = 10000

  simVals = matrix(rchisq(nsims * nEival, df=1), nrow=nsims, ncol=nEival)
  simVals = simVals %*% eB$values

  #pVal3 = mean(simVals < T0)
  pVal3 = mean(simVals > T0)
  #if( pVal3 < (50/nsims) )
  if( pVal3 == 0 )
    T3 = NA
  else
    T3 = qchisq(pVal3, df0, lower.tail = FALSE)

  dfTest = data.frame(kappa   = c(1, kappa1, kappa2, 1),
                      chiStat = c(T0, T1, T2, T3),
                      df      = c(df0, rank, d2, df0),
                      pValue  = c(pVal0, pVal1, pVal2, pVal3))
  
  row.names(dfTest) = c("Un-adjusted", "Mean scaled", "Mean-Variance scaled", "Simulated")

  .printInfoLine("Testing model fit", "Done", 50, 2)

  return(list(parEst  = dfPar,
              fitStat = dfTest))
}
