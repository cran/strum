#include "utils.h"

void process_rk(const NumericVector& yk, const NumericMatrix& xk,
                const mat& ba, mat& rka, int t)
{
  mat yka(const_cast<NumericVector&>(yk).begin(), xk.rows(), 1, false);
  mat xka(const_cast<NumericMatrix&>(xk).begin(), xk.rows(), ba.n_rows, false);

  rka = yka - (xka * ba.col(t));

  return;
}

void process_rk2(const List& yL,  const List& xL, const mat& ba,
                 const List& y2L, const List& x2L,
                 int k, int& nk, int& n2k, mat& xka, mat& rka)
{
  NumericVector yk(as<SEXP>(yL[k]));
  NumericMatrix xk(as<SEXP>(xL[k]));

  nk = xk.nrow();

  process_rk(yk, xk, ba, rka, 0);

  mat xka_(xk.begin(), xk.rows(), ba.n_rows, true);

  if( y2L.size() )
  {
    NumericVector y2k(as<SEXP>(y2L[k]));
    NumericMatrix x2k(as<SEXP>(x2L[k]));

    n2k = x2k.nrow();

    mat r2ka;

    process_rk(y2k, x2k, ba, r2ka, 1);

    rka = join_cols(rka, r2ka);

    mat x2ka(x2k.begin(), x2k.rows(), ba.n_rows, false);
    mat x12k, x21k;
    x12k.zeros(nk, ba.n_rows);
    x21k.zeros(n2k, ba.n_rows);

    xka_ = join_cols(join_rows(xka_, x12k), join_rows(x21k, x2ka));
  }

  xka = xka_;

  return;
}

void process_vk(const List& vck, const List& sigL, mat& vka,
                int n1k, int n2k, int t1, int t2)
{
  vka.zeros(n1k, n2k);

  for( int v = 0; v < sigL.size(); ++v )
  {
    NumericMatrix vck_v(as<SEXP>(vck[v]));
    NumericMatrix sig_v(as<SEXP>(sigL[v]));

    mat vck_va(vck_v.begin(), n1k, n2k, false);
    mat vcka = sig_v(t1, t2) * vck_va;

    vka += vcka;
  }

  return;
}

void process_vk2(const List& vcL,  const List& sigL,
                 const List& vc2L, const List& vc12L,
                 int k, int nk, int n2k, mat& vka)
{
  List vck(as<SEXP>(vcL[k]));

  process_vk(vck, sigL, vka, nk, nk, 0, 0);

  if( vc2L.size() )
  {
    List vc2k(as<SEXP>(vc2L[k]));
    List vc12k(as<SEXP>(vc12L[k]));

    mat v2ka;
    mat v12ka;

    process_vk(vc2k, sigL, v2ka, n2k, n2k, 1, 1);
    process_vk(vc12k, sigL, v12ka, nk, n2k, 0, 1);

    vka = join_cols(join_rows(vka, v12ka), join_rows(v12ka.t(), v2ka));
  }

  return;
}

void process_dvk2(const List& vcL, const List& vc2L, const List& vc12L, 
                  int k, int v, int nk, int n2k, int dv_type, mat& dvk)
{
  List vck(as<SEXP>(vcL[k]));

  NumericMatrix vck_v(as<SEXP>(vck[v]));

  mat dvk_v(vck_v.begin(), nk, nk, true);

  if( vc12L.size() )
  {
    mat::fixed<2, 2> sig;

    if( dv_type == 0 )
    {
      sig(0,0) = 1.0;
      sig(0,1) = sig(1,0) = sig(1,1) = 0.0;
    }
    else if( dv_type == 1 )
    {
      sig(0,0) = sig(0,1) = sig(1,0) = 0.0;
      sig(1,1) = 1.0;
    }
    else
    {
      sig(0,0) = sig(1,1) = 0.0;
      sig(0,1) = sig(1,0) = 1.0;
    }

    List vc2k(as<SEXP>(vc2L[k]));
    List vc12k(as<SEXP>(vc12L[k]));

    NumericMatrix vc2k_v(as<SEXP>(vc2k[v]));
    NumericMatrix vc12k_v(as<SEXP>(vc12k[v]));

    mat dv2k_v(vc2k_v.begin(), n2k, n2k, true);
    mat dv12k_v(vc12k_v.begin(), nk, n2k, true);

    dvk_v   = sig(0, 0) * dvk_v;
    dv2k_v  = sig(1, 1) * dv2k_v;
    dv12k_v = sig(0, 1) * dv12k_v;

    dvk_v = join_cols(join_rows(dvk_v, dv12k_v), join_rows(dv12k_v.t(), dv2k_v));
  }

  dvk = dvk_v;

  return;
}

mat get_vkidvk(const List& vcL, const List& vc2L, const List& vc12L,
                     const mat& vka,
                     int k, int v, int nk, int n2k, int dv_type)
{
  mat dvk;
  process_dvk2(vcL, vc2L, vc12L, k, v, nk, n2k, dv_type, dvk);

  mat vkidvk = solve(vka, dvk);

  return vkidvk;
}

RcppExport SEXP computeLL(SEXP ys,  SEXP xs,  SEXP vcs,
                          SEXP y2s, SEXP x2s, SEXP vc2s,
                          SEXP vc12s, SEXP betas, SEXP sigmas)
{
  double logLikelihood = 0.0;
  size_t N = 0;

  try
  {
    List yL(ys);
    List xL(xs);
    List vcL(vcs);

    List y2L(y2s);
    List x2L(x2s);
    List vc2L(vc2s);
    List vc12L(vc12s);

    List sigL(sigmas);
    mat ba = as<mat>(betas);

    size_t numBlock = yL.size();

    for( size_t k = 0; k < numBlock; ++k )
    {
      mat xka;
      mat rka;
      mat vka;

      int nk = 0, n2k = 0;

      process_rk2(yL, xL, ba, y2L, x2L, k, nk, n2k, xka, rka);

      process_vk2(vcL, sigL, vc2L, vc12L, k, nk, n2k, vka);

      mat q = rka.t() * solve(vka, rka);

      double logD, sign;
      log_det(logD, sign, vka);

      logLikelihood -= logD + q(0);
      N += (nk + n2k);
    }
  }
  catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  }
  catch( ... )
  {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  logLikelihood = .5 * logLikelihood - (N/2.0)*log(2.*M_PI);

  return wrap(logLikelihood);
}

RcppExport SEXP computeDeriv(SEXP ys,  SEXP xs,  SEXP vcs,
                             SEXP y2s, SEXP x2s, SEXP vc2s,
                             SEXP vc12s, SEXP betas, SEXP sigmas)
{
  mat dparMat;
  mat d2parMat;
  mat dpar1dsMat;
  mat dpar2dsMat;

  try
  {
    List yL(ys);
    List xL(xs);
    List vcL(vcs);

    List y2L(y2s);
    List x2L(x2s);
    List vc2L(vc2s);
    List vc12L(vc12s);

    List sigL(sigmas);
    mat ba = as<mat>(betas);

    size_t numVC    = sigL.size();
    size_t numCov   = ba.n_rows;
    size_t numTrt   = ba.n_cols;
    size_t numBlock = yL.size();

    mat dbetaMat(numCov, numBlock);
    mat dsigmaMat(numVC, numBlock);

    mat dbdbMat;
    mat dsdsMat;
    mat dbdsMat;
    mat ds1dsMat;
    mat ds2dsMat;

    dbdbMat.zeros(numCov, numCov);
    dbdsMat.zeros(numCov*numTrt, numVC);
    dsdsMat.zeros(numVC, numVC);

    ds1dsMat.zeros(numVC, numVC);
    ds2dsMat.zeros(numVC, numVC);

    for( size_t k = 0; k < numBlock; ++k )
    {
      mat xka;
      mat rka;
      mat vka;

      int nk = 0, n2k = 0;

      process_rk2(yL, xL, ba, y2L, x2L, k, nk, n2k, xka, rka);
      process_vk2(vcL, sigL, vc2L, vc12L, k, nk, n2k, vka);

      mat z = solve(vka, rka);

      // Univariate case only
      if( !y2L.size() )
      {
        mat dbeta = xka.t() * z;

        dbetaMat.col(k) = dbeta;      

        dbdbMat += -(xka.t() * solve(vka, xka));
      }

      for( size_t v = 0; v < numVC; ++v )
      {
        mat vkidvk = get_vkidvk(vcL, vc2L, vc12L, vka, k, v, nk, n2k, 2);

        mat dsigma = -0.5 * trace(vkidvk) + 0.5 * rka.t() * vkidvk * z;
        mat dbds   = -xka.t() * vkidvk * z;

        //double trvkidv = trace(vkidvk * vkidvk);
        double trvkidv = dot(vkidvk, vkidvk.t());

        mat dsds = 0.5 * trvkidv - rka.t() * vkidvk * vkidvk * z;

        dsigmaMat(v, k) = dsigma(0);
        dbdsMat.col(v) += dbds;
        dsdsMat(v, v)  += dsds(0);

        for( size_t v2 = v+1; v2 < numVC; ++v2 )
        {
          mat vkidv2k = get_vkidvk(vcL, vc2L, vc12L, vka, k, v2, nk, n2k, 2);

          //double trvkidv2 = trace(vkidvk * vkidv2k);
          double trvkidv2 = dot(vkidvk, vkidv2k.t());

          //dsds = 0.5 * trvkidv2 - rka.t() * vkidvk * vkidv2k * z;
          dsds = 0.5 * trvkidv2 - 0.5 * rka.t() * (vkidvk * vkidv2k + vkidv2k * vkidvk) * z;

          dsdsMat(v, v2) += dsds(0);
          dsdsMat(v2, v) += dsds(0);
        }
      }

      // Bivariate case only
      if( y2L.size() )
      {
        for( size_t v = 0; v < numVC; ++v )
        {
          mat vkidv1k = get_vkidvk(vcL, vc2L, vc12L, vka, k, v, nk, n2k, 0);
          mat vkidv2k = get_vkidvk(vcL, vc2L, vc12L, vka, k, v, nk, n2k, 1);

          for( size_t v2 = 0; v2 < numVC; ++v2 )
          {
            mat vkidv12k = get_vkidvk(vcL, vc2L, vc12L, vka, k, v2, nk, n2k, 2);

            //double trvkidv112 = trace(vkidv1k * vkidv12k);
            //double trvkidv212 = trace(vkidv2k * vkidv12k);
            double trvkidv112 = dot(vkidv1k, vkidv12k.t());
            double trvkidv212 = dot(vkidv2k, vkidv12k.t());

            //mat ds1ds = 0.5 * trvkidv112 - rka.t() * vkidv1k * vkidv12k * z;
            //mat ds2ds = 0.5 * trvkidv212 - rka.t() * vkidv2k * vkidv12k * z;
            mat ds1ds = 0.5 * trvkidv112 - 0.5 * rka.t() * (vkidv1k * vkidv12k + vkidv12k * vkidv1k) * z;
            mat ds2ds = 0.5 * trvkidv212 - 0.5 * rka.t() * (vkidv2k * vkidv12k + vkidv12k * vkidv2k) * z;

            ds1dsMat(v, v2) += ds1ds(0);
            ds2dsMat(v, v2) += ds2ds(0);
          }
        }
      }
    }

    if( y2L.size() )
    {
      dparMat  = dsigmaMat;
      d2parMat = dsdsMat;

      dpar1dsMat = join_cols(dbdsMat.rows(0, numCov-1), ds1dsMat);
      dpar2dsMat = join_cols(dbdsMat.rows(numCov, numCov*numTrt-1), ds2dsMat);
    }
    else
    {
      dparMat  = join_cols(dbetaMat, dsigmaMat);
      d2parMat = join_cols(join_rows(dbdbMat, dbdsMat), join_rows(dbdsMat.t(), dsdsMat));
    }
  }
  catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  }
  catch( ... )
  {
    ::Rf_error( "c++ exception (unknown reason)" );
  }

  return List::create(Named("dPar")  = dparMat,
                      Named("d2Par") = d2parMat,
                      Named("dP1dS") = dpar1dsMat,
                      Named("dP2dS") = dpar2dsMat);
}
