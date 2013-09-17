
#include <appl_grid/appl_grid.h>

extern "C" void evolvepdf_(const double&, const double&, double*);
extern "C" double alphaspdf_(const double& Q);


std::vector<double> vconvolute(appl::grid* g) {
  return g->vconvolute(evolvepdf_, alphaspdf_);
}

TH1D* hconvolute(appl::grid* g) {
  return g->convolute(evolvepdf_, alphaspdf_);
}
