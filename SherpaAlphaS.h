
#ifndef MODEL_Main_Running_AlphaS_H
#define MODEL_Main_Running_AlphaS_H

#include <algorithm>
#include <vector>

namespace SHERPA
{

//! Contains data for alpha_S running up to 3rd order
struct AsDataSet {
  double low_scale, high_scale;
  double as_low, as_high;
  int nf;
  double lambda2;
  double b[4];
  double beta0;
};
//! coefficients of the Beta functions for QCD
// cf. T. van Ritbergen, J. Vermaseren, S. Larin PLB 400 (1997) 379
// up to factor 4^(i+1)

class One_Running_AlphaS
{
  protected:

    int         m_order;
    int         m_nth, m_mzset;
    double      m_CF, m_CA;
    double      m_as_MZ, m_m2_MZ;
    double      m_fac;
    double      m_cutq2;

    AsDataSet* p_thresh;
    AsDataSet* p_shower;

    double Beta0(const int);
    double Beta1(const int);
    double Beta2(const int);
    double Beta3(const int);
    double Lambda2(const int);
    double ZetaOS2(const double, const double, const double, const int);
    double InvZetaOS2(const double, const double, const double, const int);
    double AlphaSLam(const double, const int);
    void   ContinueAlphaS(int &);

  public:

    // constructor
    One_Running_AlphaS(const int order, const double as_MZ, const double m2_MZ,
                       const std::vector<double>& qmasses);

    // destructor
    ~One_Running_AlphaS();

    // member functions
    double operator()(double q2);
    void   FixShowerLambda2(const double mu2, const double asmu,
                            const int nf=5, const int order=0);
    double AlphaS(const double q2, bool shower=false);

    // inline functions
    inline int    Order()                {
      return m_order;
    }
    inline double Beta0(const double q2) {
      return Beta0(Nf(q2));
    }
    inline double Beta1(const double q2) {
      return Beta1(Nf(q2));
    }
    inline double Beta2(const double q2) {
      return Beta2(Nf(q2));
    }
    inline double Beta3(const double q2) {
      return Beta3(Nf(q2));
    }
    inline double AsMZ()                 {
      return m_as_MZ;
    }
    int    Nf(const double q2);

    inline double CutQ2() {
      return m_cutq2;
    }
    inline double ScaleFactor() {
      return m_fac;
    }

    inline double BoundedAlphaS(const double &q2) {
      return (*this)(std::max(q2, m_cutq2));
    }
};


//! QCD SM Running of AlphaS
/*!
  \class Running_AlphaS
  \brief The class for the (running) strong coupling constant.

  This is an implementation of the
  <A HREF="http://131.169.91.193/cgi-bin/spiface/find/hep/www?rawcmd=f+eprint+hep-ph%2F9706430">
  strong coupling constant by K.~G.~Chetyrkin, B.~A.~Kniehl and M.~Steinhauser</A>

  The following points should be noted:
  - At the moment the running at up to three loops is tested, the four loop
    running remains to be checked.
  - Thresholds are fixed to the onshell masses of the correspondig quark.
  - In order to avoid the Landau Pole, the running is stopped as soon as
    \f$alpha_s\f$ equals one. (Below that scale \f$alpha_s\f$ is set linearly to zero.)
*/
/*!
  \var int Running_AlphaS::m_order
  The order at which running \f$alpha_s\f$ is evaluated, e.g. m_order = 0 refers to one-loop,
  m_order = 3 then signifies four-loop running.
*/
/*!
  \var int Running_AlphaS::m_nth
  The number of thresholds, i.e. of quarks. Later on, other particles should be implmented as well,
  at least at one-loop order. The only exception is the Z boson which is inserted as well -
  this is due to the fact that in SHERPA \f$\alpha_s\f$ at the Z-pole is used as defining
  parameter.
*/
/*!
  \var int Running_AlphaS::m_mzset
  The place of the Z-pole in the list of thresholds.
*/
/*!
  \var double Running_AlphaS::m_m2_MZ
  The mass squared of the Z-boson.
*/
/*!
  \var double Running_AlphaS::m_as_MZ
  \f$\alpha_s\f$ at the Z-pole.
*/
/*!
  \var double Running_AlphaS::Running_AlphaS::m_CF
  The Casimir operator of QCD in the fundamental representation, \f$C_F=4/3\f$.
*/
/*!
  \var double Running_AlphaS::Running_AlphaS::m_CA
  The Casimir operator of QCD in the adjoint representation, \f$C_A=3\f$.
*/
/*!
  \var AsDataSet * Running_AlphaS::p_thresh
  Values of $\alpha_s$ and the \f[\beta_i\f] at the various thresholds.
*/
/*!
  \fn double Running_AlphaS::Beta0(const double)
  Returns \f[\beta_0\f] at a given scale.
*/
/*!
  \fn double Running_AlphaS::Beta0(const int)
  Returns \f[\beta_0 = (33 - 2n_f)/12\f] for a given number of active flavours.
*/
/*!
  \fn double Running_AlphaS::Beta1(const int)
  Returns \f[\beta_1 = (306 - 38n_f)/48\f] for a given number of active flavours.
*/
/*!
  \fn double Running_AlphaS::Beta2(const int)
  Returns \f[\beta_2 = (2857/2 - 5033/18n_f + 325/54n_f^2)/64\f] for a given number of active flavours.
*/
/*!
  \fn double Running_AlphaS::Beta3(const int)
  Returns \f[\beta_3 = \left[ (149753/6 + 3564\zeta(3)) -
                      (1078361/162 + 6508/27\zeta(3))n_f +
                      (50065/162 + 6472/81\zeta(3))n_f^2 +
                      (1093/729)n_f^3\right]/256\f] for a given number of active flavours.
*/
/*!
  \fn double Running_AlphaS::Lambda2(const int)
  Fills in the b-coefficients, \f[b_{1,2,3} = \beta_{1,2,3}/\beta_0\,,\f]
  and caluclates a suitable \f$\Lambda_{\rm QCD}^2\f$ which is then returned.
*/
/*!
  \fn Running_AlphaS::Running_AlphaS(double,double)
  The constructor of Running_AlphaS fills in and sorts the thresholds according to their mass.
  \todo Documentation for this method.
*/

/*!
  \fn double Running_AlphaS::operator()(double)
  Returns the runnnig \f$alpha_s\f$ at a given scale.
*/
/*!
  \fn double Running_AlphaS::AlphaS(const double)
  Returns the value of the Running_AlphaS::operator()(double).
*/
/*!
  \fn int Running_AlphaS::Nf(const double)
  Returns the number of active flavours at a given scale.
*/

}



#endif
