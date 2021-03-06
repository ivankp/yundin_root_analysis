
#include <iostream>

#include <cstdlib>
#include <cmath>
using std::abs;
using std::pow;
using std::sqrt;

#include "SherpaAlphaS.h"

using namespace SHERPA;

namespace SHERPA
{

std::ostream &operator<< (std::ostream &str, AsDataSet &set)
{
  str << "scale->[" << set.low_scale << ", " << set.high_scale << "]";
  str << " as->[" << set.as_low << ", " << set.as_high << "]";
  str << " nf->" << set.nf <<" lam2->" << set.lambda2;
  str << " bet0->" << set.beta0;
  return str;
}

}

One_Running_AlphaS::One_Running_AlphaS(const int order,
                                       const double as_MZ, const double m2_MZ,
                                       const std::vector<double>& qmasses)
  : m_order(order),
    m_as_MZ(as_MZ), m_m2_MZ(m2_MZ), m_fac(1.),
    p_shower(NULL)
{
  p_thresh = NULL;

  m_CF    = 4./3.;
  m_CA    = 3.;

  //------------------------------------------------------------
  // SM thresholds for strong interactions, i.e. QCD
  //------------------------------------------------------------
  m_nth = qmasses.size();

  p_thresh       = new AsDataSet[m_nth+1];
  double* masses = new double[m_nth];

  for (unsigned i=0; i<qmasses.size(); i++) {
    masses[i] = qmasses[i]*qmasses[i];
  }

  std::vector<double> sortmass(&masses[0], &masses[m_nth]);
  std::sort(sortmass.begin(), sortmass.end(), std::less<double>());
  for (int i=0; i<m_nth; ++i) {
    masses[i] = sortmass[i];
  }

  int j   = 0;
  m_mzset = 0;
  for (int i=0; i<m_nth; ++j) {
    if ((masses[i] > m_m2_MZ) && (!m_mzset)) {
      //insert Z boson (starting point for any evaluation)
      m_mzset               = j;
      p_thresh[j].low_scale = m_m2_MZ;
      p_thresh[j].as_low    = m_as_MZ;
      p_thresh[j].nf        = i-1;
    }
    else {
      p_thresh[j].low_scale = masses[i];
      p_thresh[j].as_low    = 0.;
      p_thresh[j].nf        = i;
      ++i;
    }
    if (j>0) {
      p_thresh[j-1].high_scale = p_thresh[j].low_scale;
      p_thresh[j-1].as_high    = p_thresh[j].as_low;
    }
  }
  if (!m_mzset) {
    j                        = m_nth;
    m_mzset                  = j;
    p_thresh[j].low_scale    = m_m2_MZ;
    p_thresh[j].as_low       = m_as_MZ;
    p_thresh[j-1].high_scale = m_m2_MZ;
    p_thresh[j-1].as_high    = m_as_MZ;
    p_thresh[j].nf           = p_thresh[j-1].nf;
  }
  p_thresh[m_nth].high_scale   = 1.e20;
  p_thresh[m_nth].as_high      = 0.;

  delete[] masses;

  for (int i=m_mzset; i<=m_nth; ++i) {
    Lambda2(i);
    p_thresh[i].as_high       = AlphaSLam(p_thresh[i].high_scale, i);
    if (i<m_nth) {
      p_thresh[i+1].as_low    = p_thresh[i].as_high *
                                InvZetaOS2(p_thresh[i].as_high, p_thresh[i].high_scale,
                                           p_thresh[i].high_scale, p_thresh[i].nf);
    }
  }
  for (int i=m_mzset-1; i>=0; --i) {
    double lam2               = Lambda2(i);
    p_thresh[i].as_low        = AlphaSLam(p_thresh[i].low_scale, i);
    if ((lam2>p_thresh[i].low_scale) || (p_thresh[i].as_low>1.)) {
      ContinueAlphaS(i);
    } else {
      if (i>0) {
        p_thresh[i-1].as_high = p_thresh[i].as_low *
                                ZetaOS2(p_thresh[i].as_low, p_thresh[i].low_scale,
                                        p_thresh[i].low_scale, p_thresh[i-1].nf);
      }
    }
  }
}

One_Running_AlphaS::~One_Running_AlphaS()
{
  if (p_thresh != 0) {
    delete[] p_thresh;
    p_thresh = NULL;
  }
}

double One_Running_AlphaS::Beta0(const int nf)
{
  return 1./4. * (11. - (2./3.)*nf);
}

double One_Running_AlphaS::Beta1(const int nf)
{
  return 1./16. * (102. - (38./3.)*nf);
}

double One_Running_AlphaS::Beta2(const int nf)
{
  return 1./64. * (2857./2. - (5033./18.)*nf + (325./54.)*nf*nf);
}

double One_Running_AlphaS::Beta3(const int nf)
{
  double zeta3 = 1.2020569031595942854;
  return 1./256. * ( (149753./6. + 3564.*zeta3) +
                     (-1078361./162. -6508./27.*zeta3)*nf +
                     (50065./162. +6472./81.*zeta3)*(nf*nf) +
                     (1093/729)*nf*nf*nf);
}

double One_Running_AlphaS::Lambda2(const int nr)
{
  double as  = p_thresh[nr].as_low;
  double mu2 = p_thresh[nr].low_scale;
  if (as == 0.) {
    as  = p_thresh[nr].as_high;
    mu2 = p_thresh[nr].high_scale;
  }

  const double a   = as/M_PI;

  int    & nf      = p_thresh[nr].nf;
  double & beta0   = p_thresh[nr].beta0;
  double * b       = p_thresh[nr].b;
  double & lambda2 = p_thresh[nr].lambda2;

  // calculate beta coefficients
  beta0 = Beta0(nf);
  b[1]  = Beta1(nf)/beta0;
  b[2]  = Beta2(nf)/beta0;
  b[3]  = Beta3(nf)/beta0;

  double betaL = 1./a;
  if (m_order >= 1) {
    betaL     += b[1]*log(a);
    if (m_order >= 2) {
      betaL   += (b[2]-b[1]*b[1])*a;
      if (m_order >= 3) {
        betaL += (b[3]/2. - b[1] * b[2] + b[1]*b[1]*b[1]/2.)*a*a;
      }
    }
  }

  lambda2         = exp(-betaL/beta0)*mu2;
  double tas1     = AlphaSLam(mu2, nr);
  double dlambda2 = 1.e-8;
  if (abs(tas1-as)/as > 1.e-11) {
    for (; (abs(tas1-as)/as > 1.e-11);) {
      lambda2     = lambda2+dlambda2;
      double tas2 = AlphaSLam(mu2, nr);
      dlambda2    = (as-tas2)/(tas2-tas1)*dlambda2;
      tas1        = tas2;
    }
  }

  return lambda2;
}

void One_Running_AlphaS::
FixShowerLambda2(const double mu2, const double asmu,
                 const int nf, const int order)
{
  if (p_shower) {
    delete p_shower;
  }
  p_shower = new AsDataSet;
  p_shower->high_scale = mu2;
  p_shower->as_high    = asmu;
  p_shower->nf         = nf;
  p_shower->beta0      = Beta0(nf);
  p_shower->b[1]       = Beta1(nf)/p_shower->beta0;
  p_shower->b[2]       = Beta2(nf)/p_shower->beta0;
  p_shower->b[3]       = Beta3(nf)/p_shower->beta0;
  p_shower->lambda2    = 0.;

  const double a = asmu/M_PI;
  double betaL = 1./a;
  if (order >= 1) {
    betaL     += p_shower->b[1]*log(a);
    if (order >= 2) {
      betaL   += (p_shower->b[2]-p_shower->b[1]*p_shower->b[1])*a;
      if (m_order >= 3) {
        betaL += (p_shower->b[3]/2. -
                  p_shower->b[1] * p_shower->b[2] +
                  p_shower->b[1]*p_shower->b[1]*p_shower->b[1]/2.)*a*a;
      }
    }
  }

  p_shower->lambda2 = exp(-betaL/p_shower->beta0)*mu2;
  double tas1       = AlphaSLam(mu2, -1), tas2;
  double dlambda2   = 1.e-8;
  if (abs(tas1-p_shower->as_high)/p_shower->as_high > 1.e-11) {
    for (; (abs(tas1-p_shower->as_high)/p_shower->as_high > 1.e-11);) {
      p_shower->lambda2 += dlambda2;
      tas2     = AlphaSLam(mu2, -1);
      dlambda2 = (p_shower->as_high-tas2)/(tas2-tas1)*dlambda2;
      tas1     = tas2;
    }
  }
}

double One_Running_AlphaS::AlphaSLam(const double Q2, const int nr)
{
  // using shorter names
  double beta0, * b, lambda2;
  if (nr > 0) {
    beta0   = p_thresh[nr].beta0;
    b       = p_thresh[nr].b;
    lambda2 = p_thresh[nr].lambda2;
  }
  else {
    beta0   = p_shower->beta0;
    b       = p_shower->b;
    lambda2 = p_shower->lambda2;
  }

  double L = log(Q2/lambda2);
  double pref = 1./(beta0*L);
  double a = pref;
  if (m_order >= 1) {
    double logL = log(L);
    pref *= 1./(beta0*L);
    a    += -pref*(b[1] * logL);
    if (m_order >= 2) {
      double log2L = logL*logL;
      pref *= 1./(beta0*L);
      a    += pref*(b[1]*b[1]*(log2L-logL-1.) + b[2]);
      if (m_order >= 3) {
        // 3rd order (four loop) to be checked.
        double log3L = logL*log2L;
        pref *= 1./(beta0*L);
        a    += pref*(b[1]*b[1]*b[1]*(-log3L+2.5*log2L+2.*logL-0.5)
                      - 3.*b[1]*b[2] + 0.5*b[3]);
      }
    }
  }
  return M_PI*a;
}

double One_Running_AlphaS::ZetaOS2(const double as, const double mass2_os,
                                   const double mu2, const int nl)
{
  double zeta2g = 1.;

  // 0th order
  if (m_order == 0) {
    return zeta2g;
  }

  // 1st order (one loop) corresponds to two loop lambda
  double L      = log(mu2/mass2_os);
  double a      = as/M_PI;
  zeta2g       += - a*1./6.*L;
  if (m_order == 1) {
    return zeta2g;
  }

  // 2nd order
  double L2     = L*L;
  double a2     = a*a;
  zeta2g       += a2 *( 1./36.*L2 - 19./24.*L -7./24.);
  if (m_order == 2) {
    return zeta2g;
  }

  // 3rd order : not yet checked ...
  double L3     = L2*L;
  double a3     = a2*a;
  double zeta2  = M_PI*M_PI/6.;
  double zeta3  = 1.2020569031595942854;
  zeta2g       += a3 * (-58933./124416. - 2./3.*zeta2*(1.+1./3.* log(2.))
                        - 80507./27648.*zeta3 - 8521./1728.*L- 131./576. * L2
                        - 1./216.*L3 +
                        nl*(2479./31104.+ zeta2/9. + 409./1728. * L ));
  return zeta2g;
}

double One_Running_AlphaS::InvZetaOS2(const double as, const double mass2_os,
                                      const double mu2, const int nl)
{
  // might be simplified considerably when using mu2==mass2
  double zeta2g  = 1.;
  // 0th order
  if (m_order == 0) {
    return zeta2g;
  }

  // 1st order (one loop) corresponds to two loop lambda
  double L      = log(mu2/mass2_os);
  double a      = as/M_PI;
  zeta2g       += + a*1./6.*L;
  if (m_order == 1) {
    return zeta2g;
  }

  // 2nd order
  double L2     = L*L;
  double a2     = a*a;
  zeta2g       += a2 *( 1./36.*L2 + 19./24.*L + 7./24.);
  if (m_order == 2) {
    return zeta2g;
  }

  // 3rd order yet to be checked...
  double L3     = L2*L;
  double a3     = a2*a;
  double zeta2  = M_PI*M_PI/6.;
  double zeta3  = 1.2020569031595942854;
  zeta2g       += a3 * (58933./124416. + 2./3.*zeta2*(1.+1./3.* log(2.))
                        + 80507./27648.*zeta3 + 8941./1728.*L + 511./576. * L2
                        + 1./216.*L3 +
                        nl*(-2479./31104.- zeta2/9. - 409./1728. * L ));
  return zeta2g;
}



void One_Running_AlphaS::ContinueAlphaS(int & nr)
{
  // shrink actual domain
  //  * to given t0        or
  //  * to alphaS=alphaCut
  double alpha_cut = 1.;
  double & beta0   = p_thresh[nr].beta0;
  double & lambda2 = p_thresh[nr].lambda2;
  double t0        = lambda2 * exp(M_PI/(alpha_cut*beta0));
  double as        = AlphaSLam(t0, nr);
  for (; abs(as-alpha_cut) > 1.e-8;) {
    double t1      = t0 + 0.00001;
    double as1     = AlphaSLam(t1, nr);
    double das     = (as -as1)/(t0-t1);
    t1             = (alpha_cut-as)/das + t0;
    t0             = t1;
    as             = AlphaSLam(t0, nr);
  }

  m_cutq2 = t0;

  // modify lower domains
  p_thresh[nr].low_scale    = t0;
  p_thresh[nr-1].high_scale = t0;
  p_thresh[nr].as_low       = as;
  p_thresh[nr-1].as_high    = as;

  for (int i = nr-1; i>=0; --i) {
    p_thresh[i].nf          = -1;  // i.e. no ordinary running !!!
    p_thresh[i].lambda2     = 0.;
    p_thresh[i].as_low      = p_thresh[i].as_high/p_thresh[i].high_scale*p_thresh[i].low_scale;
    if (i > 0) {
      p_thresh[i-1].as_high = p_thresh[i].as_low;
    }
  }
  nr =0;
}



double One_Running_AlphaS::operator()(double q2)
{
  double as;
  q2 = q2*m_fac;
  if (q2 < 0.) {
    q2=-q2;
  }
  int i = m_mzset - 1;
  if (q2 <= m_m2_MZ) {
    for (; !((p_thresh[i].low_scale<q2) && (q2<=p_thresh[i].high_scale)); --i) {
      if (i <= 0) {
        break;
      }
    }
    if (p_thresh[i].nf >= 0) {
      as = AlphaSLam(q2, i);
    } else {
      as = q2/p_thresh[i].high_scale * p_thresh[i].as_high;
    }
  }
  else {
    ++i;
    for (; !((p_thresh[i].low_scale<q2) && (q2<=p_thresh[i].high_scale)); ++i) {
      if (i >= m_nth) {
        break;
      }
    }
    as   = AlphaSLam(q2, i);
  }
  return as;
}

double  One_Running_AlphaS::AlphaS(const double q2, bool shower)
{
  if (shower && p_shower) {
    return AlphaSLam(q2, -1);
  }
  return operator()(q2);
}

int One_Running_AlphaS::Nf(const double sc)
{
  double q2 = sc*m_fac;
  for (int i=0; i<=m_nth; ++i) {
    if (q2<=p_thresh[i].high_scale && q2>p_thresh[i].low_scale) {
      return p_thresh[i].nf;
    }
  }
  return m_nth;
}
