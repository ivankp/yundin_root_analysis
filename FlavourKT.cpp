
#include "FlavourKT.h"

using namespace fastjet;

//----------------------------------------------------------------------

class FlavourKTBJ
{
  public:
    void init(const PseudoJet& jet) {
      if (jet.E() == 0.) { // beam representative
        pt2 = 0.;
        rap = 0.;
        phi = 0.;
      } else {
        pt2 = jet.pt2();
        rap = jet.rapidity();
        phi = jet.phi();
      }
      lhid = FlavourKTPlugin::getFlavour(jet);
    }

    double distance(const FlavourKTBJ* other) const {
      // check that flavours are compatible
      const int newlhid = FlavourKTPlugin::combineFlavour(lhid, other->lhid);
      if (newlhid == 999) {
        return std::numeric_limits<double>::max();
      }
      // cannot fold beams on each other
      if (pt2 == 0. and other->pt2 == 0.) {
        return std::numeric_limits<double>::max();
      }

      double dij = 0.;
      if (pt2 == 0. or other->pt2 == 0.) {
        dij = std::max(pt2, other->pt2);
      } else {
        double dphi = abs(phi - other->phi);
        if (dphi > M_PI) {
          dphi = 2.*M_PI - dphi;
        }
        const double drap = rap - other->rap;
        dij = std::min(pt2, other->pt2)*(drap*drap + dphi*dphi);
      }
      return dij;
    }

    double beam_distance() const {
      if (pt2 == 0.) {
        return 0.75*std::numeric_limits<double>::max();
      } else {
        return 0.50*std::numeric_limits<double>::max();
      }
    }

  protected:
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
    double dij = nnh.dij_min(i, j);  // i,j are return values

    std::cout << "dij = " << dij << std::endl;

    const PseudoJet& jet1 = cs.jets()[i];
    if (not isQCD(getFlavour(jet1))) {
      dij = 0;
    }
    if (j >= 0) {
      const PseudoJet& jet2 = cs.jets()[j];
      if (not isQCD(getFlavour(jet2))) {
        dij = 0;
      }
      cs.plugin_record_ij_recombination(i, j, dij, combine(jet1, jet2), k);
      nnh.merge_jets(i, j, cs.jets()[k], k);
      std::cout << getFlavour(cs.jets()[k]) << " -> " << getFlavour(jet1) << " + " << getFlavour(jet2) << std::endl;
    } else {
      std::cout << getFlavour(jet1) << " -> beam " << std::endl;
      if (dij > 1e300) {
        dij = 0.;
      }
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
  PseudoJet newjet = PseudoJet(0., 0., 0., 0.);
  if (jet1.E() != 0. and jet2.E() != 0.) {
    newjet = jet1 + jet2;
  }
  addFlavour(newjet, combineFlavour(getFlavour(jet1), getFlavour(jet2)));
  return newjet;
}

int FlavourKTPlugin::combineFlavour(int lh1, int lh2)
{
  bool q1 = isQCD(lh1);
  bool q2 = isQCD(lh2);

  if (lh2 == 21) {
    std::swap(lh1, lh2);
    std::swap(q1, q2);
  }
  if (lh1 == 21) {
    if (q2 or lh2 == 25) {
      return lh2;
    } else {
      return 999;
    }
  }

  if (q2) {
    std::swap(lh1, lh2);
    std::swap(q1, q2);
  }
  if (q1) { // quark/anti-quark
    if (lh1 == -lh2) {
      return 21;
    } else if (not q2 and lh2 != 25) {
      return lh1;
    } else {
      return 999;
    }
  }

  return 999;
}

#if defined(__MAKECINT__)
#pragma link C++ class fastjet::JetDefinition;
#pragma link C++ class FlavourKTPlugin;
#pragma link C++ class FlavourKTBJ;
#endif
