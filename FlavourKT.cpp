
#include "FlavourKT.h"

using namespace fastjet;

//----------------------------------------------------------------------

class FlavourKTBJ
{
  public:
    void init(const PseudoJet& jet) {
      pt2 = jet.pt2();
      rap = jet.rapidity();
      phi = jet.phi();
      lhid = FlavourKTPlugin::getFlavour(jet);
    }

    double distance(const FlavourKTBJ* other) const {
      if (lhid != 21 and other->lhid != 21 and lhid != -other->lhid) {
        return std::numeric_limits<double>::max();
      }
      double dphi = abs(phi - other->phi);
      if (dphi > M_PI) {
        dphi = 2.*M_PI - dphi;
      }
      const double drap = rap - other->rap;

      double dij = std::min(pt2, other->pt2)*(drap*drap + dphi*dphi);
      return dij;
    }

    double beam_distance() const {
//       return pt2;  // that would be kT
      return 0.5*std::numeric_limits<double>::max();
    }

  private:
    double pt2;
    double rap;
    double phi;
    int lhid;
};

std::string FlavourKTPlugin::description() const
{
  return std::string("FlavourKTPlugin: if you don't know what it is, don't use it.");
}

void FlavourKTPlugin::run_clustering(ClusterSequence& cs) const
{
  int njets = cs.jets().size();
  NNH<FlavourKTBJ> nnh(cs.jets());
//   double Q2 = cs.Q2();

  while (njets > 0) {
    int i, j, k;
    double dij = nnh.dij_min(i, j); // i,j are return values...

    std::cout << "dij = " << dij << std::endl;

    if (j >= 0) {
//       cs.plugin_record_ij_recombination(i, j, dij, k);
      cs.plugin_record_ij_recombination(i, j, dij, combine(cs.jets()[i], cs.jets()[j]), k);
      nnh.merge_jets(i, j, cs.jets()[k], k);
      std::cout << getFlavour(cs.jets()[k]) << " -> " << getFlavour(cs.jets()[i]) << " + " << getFlavour(cs.jets()[j]) << std::endl;
    } else {
      cs.plugin_record_iB_recombination(i, dij);
      nnh.remove_jet(i);
    }

    njets--;
  }
}

void FlavourKTPlugin::addFlavour(PseudoJet& jet, int lhid)
{
//   FlavourInfo* finfo = new FlavourInfo(lhid);
//   jet.set_user_info(finfo);
  jet.set_user_index(lhid);
}

int FlavourKTPlugin::getFlavour(const PseudoJet& jet)
{
  return jet.user_index();
}

PseudoJet FlavourKTPlugin::combine(const PseudoJet& jet1, const PseudoJet& jet2)
{
  PseudoJet newjet = jet1 + jet2;
  addFlavour(newjet, combineFlavour(getFlavour(jet1), getFlavour(jet2)));
  return newjet;
}

int FlavourKTPlugin::combineFlavour(int lh1, int lh2)
{
  if (lh1 == 21 and lh2 == 21) {
    return 21;
  }
  if (lh1 != 21 and lh2 != 21) {
    assert(lh1 == -lh2);
    return 21;
  }
  if (lh1 != 21) {
    assert(lh2 == 21);
    return lh1;
  }
  if (lh2 != 21) {
    assert(lh1 == 21);
    return lh2;
  }
  return 9999;
}

#if defined(__MAKECINT__)
#pragma link C++ class fastjet::JetDefinition;
#pragma link C++ class FlavourInfo;
#pragma link C++ class FlavourKTPlugin;
#pragma link C++ class FlavourKTBJ;
#endif
