#ifndef TopQuarkAnalysis_TopPairBSM_interface_HTTTopJetHelper_h
#define TopQuarkAnalysis_TopPairBSM_interface_HTTTopJetHelper_h


// \class HTTTopJetHelper
// 
// \short Create tag info properties for HTTTopTags that can be computed
//        "on the fly". 
// 
//
// \author Salvatore Rappoccio
// \version first version on 1-May-2011

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/HTTTopJetTagInfo.h"

class HTTTopJetHelper : public std::unary_function<reco::Jet, reco::HTTTopJetProperties> {
 public:

  HTTTopJetHelper(double TopMass, double WMass) :
  TopMass_(TopMass), WMass_(WMass)
  {}

  reco::HTTTopJetProperties operator()( reco::Jet const & ihardJet ) const;
  
 protected:
  double      TopMass_;
  double      WMass_;

};


#endif
