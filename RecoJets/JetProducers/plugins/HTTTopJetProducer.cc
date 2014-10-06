#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/JetReco/interface/HTTTopJetTagInfo.h"
#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "HTTTopJetProducer.h"


using namespace edm;
using namespace cms;
using namespace reco;
using namespace std;

HTTTopJetProducer::HTTTopJetProducer(edm::ParameterSet const& conf):
       FastjetJetProducer( conf ),
       ptMin_(conf.getParameter<double>("jetPtMin")),
       centralEtaCut_(conf.getParameter<double>("centralEtaCut")),
       verbose_(conf.getParameter<bool>("verbose"))
{

  produces<HTTTopJetTagInfoCollection>();
  fjHEPTopTagger_ = std::auto_ptr<fastjet::HEPTopTagger>(
							 new fastjet::HEPTopTagger(conf.getParameter<double>("muCut"),
										   conf.getParameter<double>("maxSubjetMass"),
										   conf.getParameter<bool>("useSubjetMass")
										   )
							 );
  fromHTTTopJetProducer_ = 1;

}
		



void HTTTopJetProducer::produce(  edm::Event & e, const edm::EventSetup & c ) 
{
  FastjetJetProducer::produce(e, c);
}

void HTTTopJetProducer::runAlgorithm( edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if ( !doAreaFastjet_ && !doRhoFastjet_) {
    fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequence( fjInputs_, *fjJetDefinition_ ) );
  } else if (voronoiRfact_ <= 0) {
    fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequenceArea( fjInputs_, *fjJetDefinition_ , *fjAreaDefinition_ ) );
  } else {
    fjClusterSeq_ = ClusterSequencePtr( new fastjet::ClusterSequenceVoronoiArea( fjInputs_, *fjJetDefinition_ , fastjet::VoronoiAreaSpec(voronoiRfact_) ) );
  }

  //Run the jet clustering
  vector<fastjet::PseudoJet> inclusiveJets = fjClusterSeq_->inclusive_jets(ptMin_);

  if ( verbose_ ) cout << "Getting central jets" << endl;
  // Find the transient central jets
  vector<fastjet::PseudoJet> centralJets;
  for (unsigned int i = 0; i < inclusiveJets.size(); i++) {
    
    if (inclusiveJets[i].perp() > ptMin_ && fabs(inclusiveJets[i].rapidity()) < centralEtaCut_) {
      centralJets.push_back(inclusiveJets[i]);
    }
  }

  fastjet::HEPTopTagger & HEPTagger = *fjHEPTopTagger_;

  vector<fastjet::PseudoJet>::iterator jetIt = centralJets.begin(), centralJetsEnd = centralJets.end();
  if ( verbose_ )cout<<"Loop over jets"<<endl;
  for ( ; jetIt != centralJetsEnd; ++jetIt ) {
    
    if (verbose_) cout << "CMS FJ jet pt: " << (*jetIt).perp() << endl;
    
    fastjet::PseudoJet taggedJet;
    taggedJet = HEPTagger.result(*jetIt);
    
    if (taggedJet != 0){
      fjJets_.push_back(taggedJet);           
    }
  }
  
}

void HTTTopJetProducer::addHTTTopJetTagInfoCollection( edm::Event& iEvent, 
						       const edm::EventSetup& iSetup,
						       edm::OrphanHandle<reco::BasicJetCollection> & oh){


  // Set up output list
  auto_ptr<HTTTopJetTagInfoCollection> tagInfos(new HTTTopJetTagInfoCollection() );

  // Loop over jets
  for (size_t ij=0; ij != fjJets_.size(); ij++){

    HTTTopJetProperties properties;
    HTTTopJetTagInfo tagInfo;

    // Black magic:
    // In the standard CA treatment the RefToBase is made from the handle directly
    // Since we only have a OrphanHandle (the JetCollection is created by this process) 
    // we have to take the detour via the Ref
    edm::Ref<reco::BasicJetCollection> ref(oh, ij);  
    edm::RefToBase<reco::Jet> rtb(ref);  
    
    properties.topMass          = 0.;
    properties.rawMass	        = 0.;
    properties.prunedMass	= 0.;
    properties.fW		= 0.;
    properties.massRatioPassed  = 0.;

    // Only needed for MultiR tagger
    properties.isMultiR	        = 0;
    properties.Rmin	        = -1.;
    properties.RminExpected     = -1.;
    
    tagInfo.insert(rtb, properties );
    tagInfos->push_back( tagInfo );
  }  

  iEvent.put( tagInfos );
  
};

 
//define this as a plug-in
DEFINE_FWK_MODULE(HTTTopJetProducer);
