#ifndef GENMATCHFITPRODUCER_H
#define GENMATCHFITPRODUCER_H

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "SRothman/DataFormats/interface/EEC.h"

#include <iostream>
#include <memory>
#include <vector>

using vecptr_d_t = std::shared_ptr<std::vector<double>>;

class GenMatchFitProducer : public edm::stream::EDProducer<> {
public:
  explicit GenMatchFitProducer(const edm::ParameterSet&);
  ~GenMatchFitProducer() override {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::InputTag jetsTag_;
  edm::EDGetTokenT<EECPartsCollection> jetsToken_;

  edm::InputTag genJetsTag_;
  edm::EDGetTokenT<EECPartsCollection> genJetsToken_;

  double dR2cut_;
  double minPartPt_;
  double partDR2cut_;

  vecptr_d_t recoCorrPT_, errPT_, errETA_, errPHI_;

  unsigned maxIter_;
  double feasCondition_;
  double startMu_;
  double startLambda_;
  double clipVal_;
};

#endif
