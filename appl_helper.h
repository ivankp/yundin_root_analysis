
#include <appl_grid/appl_grid.h>

extern "C" void evolvepdf_(const double&, const double&, double*);
extern "C" double alphaspdf_(const double& Q);


std::vector<double> vconvolute(appl::grid* g, int nloops=0,
                               double ren_scl_fac=1., double fac_scl_fac=1.,
                               double energy_fac=1.)
{
  return g->vconvolute(evolvepdf_, alphaspdf_, nloops, ren_scl_fac, fac_scl_fac, energy_fac);
}

// TH1D* hconvolute(appl::grid* g, int nloops=0,
//                  double ren_scl_fac=1., double fac_scl_fac=1.,
//                  double energy_fac=1.)
// {
//   return g->convolute(evolvepdf_, alphaspdf_, nloops, ren_scl_fac, fac_scl_fac, energy_fac);
// }

std::vector<double> pvconvolute(int subproc, appl::grid* g, int nloops=0,
                               double ren_scl_fac=1., double fac_scl_fac=1.,
                               double energy_fac=1.)
{
  return g->vconvolute_subproc(subproc, evolvepdf_, alphaspdf_, nloops, ren_scl_fac, fac_scl_fac, energy_fac);
}
