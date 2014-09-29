#include "HTTTopJetTagger.h"
#include "RecoJets/JetAlgorithms/interface/HTTTopJetHelper.h"
#include "DataFormats/JetReco/interface/HTTTopJetTagInfo.h"

using namespace std;
using namespace reco;
using namespace edm;
//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
HTTTopJetTagger::HTTTopJetTagger(const edm::ParameterSet& iConfig):
  src_(iConfig.getParameter<InputTag>("src") ),
  TopMass_(iConfig.getParameter<double>("TopMass") ),
  WMass_(iConfig.getParameter<double>("WMass") ),
  verbose_(iConfig.getParameter<bool>("verbose") )
{
  produces<HTTTopJetTagInfoCollection>();
  
  input_jet_token_ = consumes<edm::View<reco::Jet> >(src_);

}


HTTTopJetTagger::~HTTTopJetTagger()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HTTTopJetTagger::produce( edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Set up output list
  auto_ptr<HTTTopJetTagInfoCollection> tagInfos(new HTTTopJetTagInfoCollection() );

  // Get the input list of basic jets corresponding to the hard jets
  Handle<View<Jet> > pBasicJets;
  iEvent.getByToken(input_jet_token_, pBasicJets);

  // Get a convenient handle
  View<Jet> const & hardJets = *pBasicJets;

  HTTTopJetHelper helper( TopMass_, WMass_ );
   
  // Now loop over the hard jets and do kinematic cuts
  View<Jet>::const_iterator ihardJet = hardJets.begin(),
    ihardJetEnd = hardJets.end();
  size_t iihardJet = 0;
  for ( ; ihardJet != ihardJetEnd; ++ihardJet, ++iihardJet ) {

    if ( verbose_ ) cout << "Processing ihardJet with pt = " << ihardJet->pt() << endl;

    // Initialize output variables
    // Get a ref to the hard jet
    RefToBase<Jet> ref( pBasicJets, iihardJet );    
    // Get properties
    HTTTopJetProperties properties = helper( *ihardJet );
    
    HTTTopJetTagInfo tagInfo;
    tagInfo.insert( ref, properties );
    tagInfos->push_back( tagInfo );
  }// end loop over hard jets
  
  iEvent.put( tagInfos );
 
  return;   
}


// ------------ method called once each job just before starting event loop  ------------
void 
HTTTopJetTagger::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HTTTopJetTagger::endJob() {

}

//define this as a plug-in
DEFINE_FWK_MODULE(HTTTopJetTagger);
