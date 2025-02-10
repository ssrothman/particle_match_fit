#ifndef PMF_CHISQ_LOSS_FCN_H
#define PMF_CHISQ_LOSS_FCN_H

#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/SimonTools/src/util.h"

#include "Minuit2/FCNBase.h"
#include <vector>
#include <memory>

enum spatialLoss{
    TYPE1=0, //construct pT-weighted predicted pT, eta. Use in classic chisq loss
    TYPE2=1  //weight the chisq difference between recoETA, genETA by transfer matrix
};

class ChisqLossFCN: public ROOT::Minuit2::FCNBase {
  private:
    //data defining fit problem
    Eigen::VectorXd recoPT, recoETA, recoPHI;
    Eigen::VectorXd genPT, genETA, genPHI;
    Eigen::VectorXd errPT, errETA, errPHI;

    size_t NPReco, NPGen;

    Eigen::VectorXd weightedGenETA, weightedGenPHI;

    std::vector<std::pair<unsigned, unsigned>> locations;

    enum spatialLoss type;

    std::vector<double> PUpt0s, PUexps, PUpenalties;

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

    explicit ChisqLossFCN(const simon::jet& recojet,
                          const simon::jet& genjet,
                          const std::vector<std::pair<unsigned, unsigned>>& locations,
                          const enum spatialLoss type,
                          const std::vector<double>& PUpt0s,
                          const std::vector<double>& PUexps, 
                          const std::vector<double>& PUpenalties);

    double operator()(const std::vector<double>& data) const override;
    
    //error computation constant
    //should be 1.0 for our chisq likelihood 
    inline double Up() const override {return 1.0;}
};

#endif
