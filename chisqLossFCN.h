#ifndef PMF_CHISQ_LOSS_FCN_H
#define PMF_CHISQ_LOSS_FCN_H

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/util.h"

#include "Minuit2/FCNBase.h"
#include <vector>
#include <memory>
#include <armadillo>

enum spatialLoss{
    TYPE1=0, //construct pT-weighted predicted pT, eta. Use in classic chisq loss
    TYPE2=1  //weight the chisq difference between recoETA, genETA by transfer matrix
};

class ChisqLossFCN: public ROOT::Minuit2::FCNBase {
  private:
    //data defining fit problem
    const arma::vec recoPT, recoETA, recoPHI;
    const arma::vec genPT, genETA, genPHI;
    const arma::vec errPT, errETA, errPHI;

    const size_t NPReco, NPGen;

    const arma::vec weightedGenETA, weightedGenPHI;

    const std::vector<std::pair<unsigned, unsigned>> locations;

    const enum spatialLoss type;

    const std::vector<double> PUpt0s, PUexps, PUpenalties;

    std::vector<unsigned> ids;

  public:
    ChisqLossFCN() : 
        recoPT(), recoETA(), recoPHI(),
        genPT(), genETA(), genPHI(),
        errPT(), errETA(), errPHI(),
        NPReco(0), NPGen(0),
        weightedGenETA(), weightedGenPHI(),
        locations(),
        type(TYPE1),
        PUpt0s(), PUexps(), PUpenalties(),
        ids() {}

    explicit ChisqLossFCN(const jet& recojet,
                          const jet& genjet,
                          const std::vector<std::pair<unsigned, unsigned>>& locations,
                          const enum spatialLoss type,
                          const std::vector<double>& PUpt0s,
                          const std::vector<double>& PUexps, 
                          const std::vector<double>& PUpenalties):
      recoPT(recojet.ptvec()), 
      recoETA(recojet.etavec()), 
      recoPHI(recojet.phivec()),
      genPT(genjet.ptvec()), 
      genETA(genjet.etavec()), 
      genPHI(genjet.phivec()),
      errPT(recojet.dptvec()), 
      errETA(recojet.detavec()), 
      errPHI(recojet.dphivec()),
      NPReco(recoPT.size()), NPGen(genPT.size()),
      weightedGenETA(genPT % genETA), 
      weightedGenPHI(genPT % genPHI),
      locations(locations),
      type(type),
      PUpt0s(PUpt0s),
      PUexps(PUexps), PUpenalties(PUpenalties),
      ids(){

          for(const auto& part : recojet.particles){
            if(part.pdgid == 22){
                ids.emplace_back(0);
            } else if(part.pdgid == 130){
                ids.emplace_back(1);
            } else if(part.charge != 0){
                ids.emplace_back(2);
            } else {
                throw std::runtime_error("Unrecognized particle type in recojet");
            }
          }
      }

    double operator()(const std::vector<double>& data) const override;
    
    //error computation constant
    //should be 1.0 for our chisq likelihood 
    inline double Up() const override {return 1.0;}
};

#endif
