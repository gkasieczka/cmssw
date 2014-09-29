#include "RecoJets/JetAlgorithms/interface/HTTTopJetHelper.h"


struct GreaterByPtCandPtr {
  bool operator()( const edm::Ptr<reco::Candidate> & t1, const edm::Ptr<reco::Candidate> & t2 ) const {
    return t1->pt() > t2->pt();
  }
};


reco::HTTTopJetProperties HTTTopJetHelper::operator()( reco::Jet const & ihardJet ) const {
  
  reco::HTTTopJetProperties properties;

  // // Get subjets  
  // reco::Jet::Constituents subjets = ihardJet.getJetConstituents();
  // properties.nSubJets = subjets.size();  // number of subjets
  // properties.topMass = ihardJet.mass();      // jet mass
  // properties.wMass = 99999.;                  // best W mass
  // properties.minMass = 999999.;            // minimum mass pairing
  
  return properties;
}
