
#ifndef FLAVOURKT_H
#define FLAVOURKT_H

#ifndef __CINT__
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/NNH.hh>
#else
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/NNH.hh"
#endif /* __CINT __ */

#include "fastjet/ClusterSequence.hh"
#include <algorithm>
#include <limits>

// class FlavourInfo : public fastjet::PseudoJet::UserInfoBase
// {
//   public:
//     FlavourInfo(int lhid_)
//       : lhid(lhid_)
//     {}
//
//     const int lhid;
// };

class FlavourKTPlugin : public fastjet::JetDefinition::Plugin
{
  friend class FlavourKTBJ;
  public:
    FlavourKTPlugin(double radius_=1.)
      : radius(radius_)
    {}

    virtual std::string description() const;
    virtual void run_clustering(fastjet::ClusterSequence& cs) const;

    virtual double R() const {
      return radius;
    }

    static int getFlavour(const fastjet::PseudoJet& jet);
    static void addFlavour(fastjet::PseudoJet& jet, int lhid);


  protected:
    static fastjet::PseudoJet combine(const fastjet::PseudoJet& jet1,
                                      const fastjet::PseudoJet& jet2);
    static int combineFlavour(int lh1, int lh2);
    static bool isQCD(int id) { return id == 21 or abs(id) <= 6; }

    double radius;
};

#endif /*FLAVOURKT_H*/
