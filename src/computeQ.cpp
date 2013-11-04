#include "utils.h"

mat make_delta(const mat& C, const List& V)
{
  mat deltaVec;

  size_t numVC  = V.size();
  size_t numCov = C.n_cols;
  size_t numTrt = C.n_rows;

  size_t numPar = (numCov+numVC)*numTrt + numTrt*(numTrt-1)/2*numVC;
  deltaVec.zeros(numPar, 1);

  size_t i = 0;

  for( size_t t1 = 0; t1 < numTrt; ++t1 )
  {
    for( size_t p = 0; p < numCov; ++p )
      deltaVec(i++, 0) = C(t1, p);

    for( size_t j = 0; j < numVC; ++j )
    {
      deltaVec(i++, 0) = as<mat>(V[j])(t1, t1);
    }
    for( int t2 = t1-1; t2 >= 0; --t2 )
    {
      for( size_t j = 0; j < numVC; ++j )
      {
        deltaVec(i++, 0) = as<mat>(V[j])(t1, t2);
      }
    }
  }

  return deltaVec;
}

RcppExport SEXP getDelta(SEXP C_, SEXP V_)
{
  mat deltaVec;

  try
  {
    mat  C = as<mat>(C_);
    List V(V_);

    deltaVec = make_delta(C, V);
  }
  catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  }
  catch( ... )
  {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return wrap(deltaVec);
}

RcppExport SEXP computeDelta(SEXP B_, SEXP L_, SEXP Gs_, SEXP Gm_, SEXP Z_, SEXP E_)
{
  mat B  = as<mat>(B_);
  mat L  = as<mat>(L_);
  mat Gs = as<mat>(Gs_);
  mat Gm = as<mat>(Gm_);
  List Z(Z_);
  List E(E_);

  mat I = eye<mat>(B.n_rows, B.n_rows);
  mat H;
  try
  {
    H = inv(I - B);
  }
  catch( std::runtime_error )
  {
    H = pinv(I - B);
  }
  catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  }
  catch( ... )
  {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  mat LH = L * H;
  mat C = Gm + LH * Gs;

  List V((size_t)Z.size());

  for( size_t j = 0; j < (size_t)Z.size(); ++j )
  {
    mat Zj = as<mat>(Z[j]);
    mat Ej = as<mat>(E[j]);

    mat Vj = LH * Zj * LH.t() + Ej;
    V[j] = Vj;
  }

  mat deltaVec = make_delta(C, V);

  return wrap(deltaVec);
}

mat process_Amat(SEXP B_, SEXP L_, SEXP Gs_, SEXP Gm_, List Zs, List Es, const mat& intercept)
{
  size_t lengthV = Es.size() - 1;

  mat As = zeros<mat>(intercept.n_rows, lengthV);

  for( size_t i = 0; i < lengthV; ++i )
  {
    mat thetaNVei = as<mat>(computeDelta(B_, L_, Gs_, Gm_, as<SEXP>(Zs[i+1]), as<SEXP>(Es[i+1])));

    As.col(i) = (thetaNVei - intercept);
  }

  return(As);
}

mat process_thetaV(const mat& delta, const mat& w, const mat& intercept, const mat& As)
{
  mat thetaV;

  mat AWA = As.t() * diagmat(w) * As;
  mat AWd = As.t() * diagmat(w) * (delta - intercept);

  try
  {
    thetaV = solve(AWA, AWd);
  }
  catch( std::runtime_error )
  {
    mat AWAi;

    try
    {
      AWAi = inv(AWA);
    }
    catch( std::runtime_error )
    {
      AWAi = pinv(AWA);
    }
    catch( std::exception& __ex__ )
    {
      forward_exception_to_r( __ex__ );
    }
    catch( ... )
    {
      ::Rf_error( "c++ exception (unknown reason)" );
    }

    thetaV = AWAi * AWd;
  }
  catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  }
  catch( ... )
  {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return(thetaV);
}

RcppExport SEXP computeThetaV(SEXP B_,  SEXP L_,  SEXP Gs_,    SEXP Gm_,
                              SEXP Zs_, SEXP Es_, SEXP delta_, SEXP w_)
{
  mat thetaV;

  try
  {
    List Zs(Zs_);
    List Es(Es_);

    mat delta = as<mat>(delta_);
    mat w     = as<mat>(w_);

    mat intercept = as<mat>(computeDelta(B_, L_, Gs_, Gm_, as<SEXP>(Zs[0]), as<SEXP>(Es[0])));
    mat aMat      = process_Amat(B_, L_, Gs_, Gm_, Zs, Es, intercept);

    thetaV = process_thetaV(delta, w, intercept, aMat);
  }
  catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  }
  catch( ... )
  {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return wrap(thetaV);
}

RcppExport SEXP computeQstar(SEXP B_,  SEXP L_,  SEXP Gs_,    SEXP Gm_,
                             SEXP Zs_, SEXP Es_, SEXP delta_, SEXP w_)
{
  double qStar = 0.0;

  try
  {
    List Zs(Zs_);
    List Es(Es_);

    mat delta = as<mat>(delta_);
    mat w     = as<mat>(w_);

    mat intercept = as<mat>(computeDelta(B_, L_, Gs_, Gm_, as<SEXP>(Zs[0]), as<SEXP>(Es[0])));
    mat aMat      = process_Amat(B_, L_, Gs_, Gm_, Zs, Es, intercept);

    mat thetaV = process_thetaV(delta, w, intercept, aMat);
    mat diff   = delta - (intercept + (aMat * thetaV));

    mat qStarM = sum(diff % diff % w);
    qStar = qStarM(0,0);
  }
  catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  }
  catch( ... )
  {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return wrap(qStar);
}
