#include "DataFormats/Math/interface/deltaR.h"
#include "SRothman/DataFormats/interface/EMDFlow.h"

#include "SRothman/armadillo/include/armadillo"


#include "SRothman/GenMatchFit/plugins/GenMatchFitProducer.h"
#include "SRothman/GenMatchFit/plugins/GenMatchFCN.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"

#define EPSILON 1e-15

using namespace ROOT::Minuit2;

// requires at least C++11
//from https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf/49812018#49812018
static const std::string vformat(const char * const zcFormat, ...) {

  // initialize use of the variable argument array
  va_list vaArgs;
  va_start(vaArgs, zcFormat);

  // reliably acquire the size
  // from a copy of the variable argument array
  // and a functionally reliable call to mock the formatting
  va_list vaArgsCopy;
  va_copy(vaArgsCopy, vaArgs);
  const int iLen = std::vsnprintf(NULL, 0, zcFormat, vaArgsCopy);
  va_end(vaArgsCopy);

  // return a formatted string without risking memory mismanagement
  // and without assuming any compiler or platform specific behavior
  std::vector<char> zc(iLen + 1);
  std::vsnprintf(zc.data(), zc.size(), zcFormat, vaArgs);
  va_end(vaArgs);
  return std::string(zc.data(), iLen); 
}

GenMatchFitProducer::GenMatchFitProducer(const edm::ParameterSet& conf)
    : 
      jetsTag_(conf.getParameter<edm::InputTag>("jets")),
      jetsToken_(consumes<EECPartsCollection>(jetsTag_)),
      genJetsTag_(conf.getParameter<edm::InputTag>("genJets")),
      genJetsToken_(consumes<EECPartsCollection>(genJetsTag_)),
      dR2cut_(conf.getParameter<double>("dR2cut")),
      minPartPt_(conf.getParameter<double>("minPartPt")),
      partDR2cut_(conf.getParameter<double>("partDR2cut")),
      recoCorrPT_(std::make_shared<std::vector<double>>()),
      errPT_(std::make_shared<std::vector<double>>()),
      errETA_(std::make_shared<std::vector<double>>()),
      errPHI_(std::make_shared<std::vector<double>>()),
      maxIter_(conf.getParameter<unsigned>("maxIter")),
      feasCondition_(conf.getParameter<double>("feasCondition")),
      startMu_(conf.getParameter<double>("startMu")),
      startLambda_(conf.getParameter<double>("startLambda")),
      clipVal_(conf.getParameter<double>("clipVal"))
{
        produces<EMDFlowCollection>();
}

void GenMatchFitProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets");
  desc.add<edm::InputTag>("genJets");
  desc.add<double>("dR2cut");
  desc.add<double>("minPartPt");
  desc.add<double>("partDR2cut");
  desc.add<unsigned>("maxIter");
  desc.add<double>("feasCondition");
  desc.add<double>("startMu");
  desc.add<double>("startLambda");
  desc.add<double>("clipVal");
  descriptions.addWithDefaultLabel(desc);
}

void GenMatchFitProducer::produce(edm::Event& evt, const edm::EventSetup& setup) {
  edm::Handle<EECPartsCollection> jets;
  evt.getByToken(jetsToken_, jets);

  edm::Handle<EECPartsCollection> genJets;
  evt.getByToken(genJetsToken_, genJets);

  auto result = std::make_unique<EMDFlowCollection>();
  
  std::vector<unsigned> taken;

  for(unsigned iReco=0; iReco<jets->size(); ++iReco){//for each reco jet
    const auto& recoJet = jets->at(iReco);

    int bestGen=-1;
    double bestDR2=100000;

    for(unsigned iGen=0; iGen<genJets->size(); ++iGen){//for each gen jet
      for(auto take : taken){//check if already matched
        if(take==iGen){
          continue;
        }
      }//end check if already matched

      const auto& genJet = genJets->at(iGen);

      double dR2 = reco::deltaR2(recoJet.jetEta, recoJet.jetPhi, 
                                 genJet.jetEta, genJet.jetPhi);
      if(dR2>dR2cut_ || dR2>bestDR2){//update best match
        continue;
      } else {
        bestGen = iGen;
        bestDR2 = dR2;
      }
    }//end for each gen jet

    if(bestGen>=0){//if matched
      const auto& genJet = genJets->at(bestGen);

      recoCorrPT_->clear();
      errPT_->clear();
      errETA_->clear();
      errPHI_->clear();
      double Efactor = genJet.jetPt/recoJet.rawPt; //ideal JEC factor
      for(size_t iPart=0; iPart<recoJet.partPt->size(); ++iPart){
        double newPt = recoJet.partPt->at(iPart)*Efactor;
        recoCorrPT_->push_back(newPt); 
        errPT_->push_back(0.02 * newPt);
        errETA_->push_back(0.05);
        errPHI_->push_back(0.05);
      }
  
      unsigned NPReco = recoJet.partPt->size();
      unsigned NPGen = genJet.partPt->size();

      GenMatchFCN FCN(recoCorrPT_, recoJet.partEta, recoJet.partPhi,
                      genJet.partPt, genJet.partEta, genJet.partPhi,
                      errPT_, errETA_, errPHI_, 
                      3.0, 1.0,
                      startMu_, startLambda_);

      MnUserParameters upar;
      for(size_t iPGen=0; iPGen<NPGen; ++iPGen){
        std::string name = vformat("lambda%lu", iPGen);
        upar.Add(name, 1.0, 1.0);
        upar.Fix(name);
      }
      upar.Add("mu", 1.0, 1.0);

      for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
        for(size_t iPGen=0; iPGen<NPGen; ++iPGen){

          double dR2 = reco::deltaR2(recoJet.partEta->at(iPReco), 
                                     recoJet.partPhi->at(iPReco),
                                     genJet.partEta->at(iPGen), 
                                     genJet.partPhi->at(iPGen));

          std::string name = vformat("%lux%lu", iPReco, iPGen);
          if(dR2 < partDR2cut_){
            upar.Add(name, 0.5, 0.5); 
            upar.SetLowerLimit(name, 0.0);
          } else{
            upar.Add(name, 0., 0.);
            upar.Fix(name);
          }
        }
      }

      size_t iIter=0;
      bool converged=false;
      double feas=0;
      MnMigrad migrad(FCN, upar);

      while(iIter++<maxIter_ && !converged){
        FunctionMinimum min = migrad();

        auto state = min.UserState();

        feas = FCN.getFeas(state.Parameters().Params());
        //FCN.updateMu(2.0); 

        if(feas < feasCondition_){
          converged=true;
        } else{
        }

        migrad.SetValue("mu", migrad.Value("mu") * 3);

        for(size_t iPGen=0; iPGen<NPGen; ++iPGen){
          double C=0;
          for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
            C += migrad.Value(FCN.idx(iPReco, iPGen));
          }
          if (C!=0){
            for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
              size_t idx = FCN.idx(iPReco, iPGen);
              migrad.SetValue(idx, migrad.Value(idx)/C);
            }
          }
        }

        for(size_t iPGen=0; iPGen<NPGen; ++iPGen){
          for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
            size_t idx = FCN.idx(iPReco, iPGen);
            double val = migrad.Value(idx);
            if(val < clipVal_){
              migrad.SetValue(idx, 0.0);
              migrad.Fix(idx);
            } else if (1-val < clipVal_){
              migrad.SetValue(idx, 1.0);
              migrad.Fix(idx);
            }
          }
        }
      }

      std::vector<double> C2(NPGen, 0);
      for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
        for(size_t iPGen=0; iPGen<NPGen; ++iPGen){
          std::string name = vformat("%lux%lu", iPReco, iPGen);
          C2.at(iPGen) += migrad.Value(name.c_str());
        }
      }

      auto EG = std::make_shared<std::vector<double>>();
      arma::vec EGvec(NPGen, arma::fill::zeros);
      auto ER = std::make_shared<std::vector<double>>();
      arma::vec ERvec(NPReco, arma::fill::zeros);
      for(unsigned iPReco=0; iPReco<NPReco; ++iPReco){
        ER->push_back(recoJet.partPt->at(iPReco) / recoJet.rawPt);
        ERvec(iPReco) = recoJet.partPt->at(iPReco) / recoJet.rawPt;
      }
      for(unsigned iPGen=0; iPGen<NPGen; ++iPGen){
        EG->push_back(genJet.partPt->at(iPGen) / genJet.jetPt);
        EGvec(iPGen) = genJet.partPt->at(iPGen) / genJet.jetPt;
      }

      auto flowmat = std::make_shared<arma::mat>(NPReco, NPGen, arma::fill::zeros);
      for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
        for(size_t iPGen=0; iPGen<NPGen; ++iPGen){
          if(C2.at(iPGen)==0){
            (*flowmat)(iPReco, iPGen) = 0;
          } else{
            std::string name = vformat("%lux%lu", iPReco, iPGen);
            (*flowmat)(iPReco, iPGen) = migrad.Value(name.c_str())/C2.at(iPGen); //* ER->at(iPReco)/EG->at(iPGen);
          }
        }
      }

      arma::vec EP = *flowmat * EGvec;
      for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
        if(EP(iPReco)==0){
          continue;
        }
        for(size_t iPGen=0; iPGen<NPGen; ++iPGen){
          (*flowmat)(iPReco, iPGen) *= ERvec(iPReco)/EP(iPReco);
        }
      }

      std::cout << "Reco particles:" << std::endl;
      for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
        printf("(%0.3f, %0.3f, %0.3f, %d, %d)\n", recoJet.partPt->at(iPReco), 
                                                  recoJet.partEta->at(iPReco),
                                                  recoJet.partPhi->at(iPReco),
                                                  recoJet.partCharge->at(iPReco),
                                                  recoJet.partPdgId->at(iPReco));
      }
      std::cout << std::endl << "Gen particles:" << std::endl;
      for(size_t iPGen=0; iPGen<NPGen; ++iPGen){
        printf("(%0.3f, %0.3f, %0.3f, %d, %d)\n", genJet.partPt->at(iPGen), 
                                                  genJet.partEta->at(iPGen),
                                                  genJet.partPhi->at(iPGen),
                                                  genJet.partCharge->at(iPGen),
                                                  genJet.partPdgId->at(iPGen));
      }

      std::cout << std::endl << "Flow matrix" << std::endl << *flowmat << std::endl;
      std::cout << "genPT" << std::endl << EGvec << std::endl;
      std::cout << "recoPT" << std::endl << ERvec << std::endl;
      std::cout << "F * EG" << std::endl << *flowmat * EGvec << std::endl;

      arma::rowvec sumG = arma::sum(*flowmat, 0);
      arma::colvec sumR = arma::sum(*flowmat, 1);

      auto matchedG = std::make_shared<std::vector<bool>>();
      auto matchedR = std::make_shared<std::vector<bool>>();

      for(size_t iPGen=0; iPGen<NPGen; ++iPGen){
        matchedG->push_back(sumG(iPGen)!=0);
      }
      for(size_t iPReco=0; iPReco<NPReco; ++iPReco){
        matchedR->push_back(sumR(iPReco!=0));
      }

      result->emplace_back(genJet.iJet, recoJet.iJet, std::move(flowmat), 
          NPGen, NPReco,
          std::move(EG), std::move(ER),
          std::move(matchedG), std::move(matchedR));
    }//end if matched
  }//end for each reco jet
  evt.put(std::move(result));
}  // end produce()

DEFINE_FWK_MODULE(GenMatchFitProducer);
