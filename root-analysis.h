#include <string>
#include <vector>

namespace RootAnalysis {
  struct NTupleEvent {
    long int id;
    int nparticle;
    float* px;
    float* py;
    float* pz;
    float* E;
    double alphas;
    int* kf;
    double weight, me_wgt;
    double x1, x2, x1p, x2p;
    int id1, id2;
    double fac_scale, ren_scale;
    int nuwgt;
    double usr_wgts[18];
    int alphaspower;
    char part[2];
    int trials;
  };

  int Init(const std::vector<std::string>& cmdline);
  int Analyse(const NTupleEvent& event);
  int Finish();
}
