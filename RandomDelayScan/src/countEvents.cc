#include "TrackerDAQAnalysis/RandomDelayScan/interface/countEvents.h"

//! ctor
countEvents::countEvents(const edm::ParameterSet& iConfig) {
  numberOfEvents_ = 0;
}

// ----------------------------------------------------------------

//! dtor
countEvents::~countEvents()
{}

// ----------------------------------------------------------------


//! loop over the reco particles and count leptons
bool countEvents::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  numberOfEvents_++;
  return true;
}

void countEvents::endJob() {
  std::cout<<"Number of events = "<<numberOfEvents_<<std::endl;
}

DEFINE_FWK_MODULE(countEvents);
