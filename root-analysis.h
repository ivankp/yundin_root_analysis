
#include <string>
#include <vector>

namespace RootAnalysis {

  struct NTupleEvent {
    int id;
    int nparticle;
    float* px;
    float* py;
    float* pz;
    float* E;
    double alphas;
    int* kf;
    double weight;
    double weight2; // unused
    double me_wgt;
    double me_wgt2; // unused
    double x1;
    double x2;
    double x1p;
    double x2p;
    int id1;
    int id2;
    double fac_scale;
    double ren_scale;
    int nuwgt;
    double usr_wgts[18];
    char part[2];
    char alphaspower;
    double trials;
  };

#if !defined(__MAKECINT__)
  bool Init(const std::vector<std::string>& cmdline, const NTupleEvent* event);
  bool Analyse(const NTupleEvent* event);
  bool Finish();
#endif
}
