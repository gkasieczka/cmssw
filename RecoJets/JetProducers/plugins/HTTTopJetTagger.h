#ifndef TopQuarkAnalysis_TopJetProducers_interface_HTTTopJetTagger_h
#define TopQuarkAnalysis_TopJetProducers_interface_HTTTopJetTagger_h

// -*- C++ -*-
//
// Package:    HTTTopJetTagger
// Class:      HTTTopJetTagger
// 
/**\class HTTTopJetTagger HTTTopJetTagger.cc TopQuarkAnalysis/TopJetProducers/src/HTTTopJetTagger.cc


*/
//
// Original Author:  "Salvatore Rappoccio"
//         Created:  Thu Jul  3 00:19:30 CDT 2008
//
//


// system include files
#include <memory>
#include <vector>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include <Math/VectorUtil.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>


//
// class decleration
//

class HTTTopJetTagger : public edm::EDProducer {
   public:
      explicit HTTTopJetTagger(const edm::ParameterSet&);
      ~HTTTopJetTagger();


   private:
      virtual void beginJob() ;
      virtual void produce( edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  edm::InputTag   src_;

  double      TopMass_;
  double      WMass_;
  bool        verbose_;

  edm::EDGetTokenT<edm::View<reco::Jet> > input_jet_token_;

};



#endif
