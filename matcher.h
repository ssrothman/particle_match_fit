#ifndef PMF_MATCHER_H
#define PMF_MATCHER_H

#include <vector>
#include <unordered_set>
#include <iostream>

#include <armadillo>

#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/util.h"

#include "chisqLossFCN.h"
#include "MatchingFilter.h"
#include "ParticleUncertainty.h"
#include "prefit.h"

using namespace ROOT::Minuit2;

class matcher{
public:
    explicit matcher(const jet& recojet,
                     const jet& genjet,
                     const std::vector<bool>& excludeGen,

                     double clipval, 

                     const enum spatialLoss& loss,
                     const enum matchFilterType& filter,
                     const enum uncertaintyType& uncertainty,
                     const std::vector<enum prefitterType>& prefitters,
                     double PUexp, double PUpenalty,

                     bool recoverLostTracks,

                     double cutoff, 

                     double softPt=0, double hardPt=1e9,

                     const std::vector<double>& EMstochastic = {}, 
                     const std::vector<double>& EMnoise = {},
                     const std::vector<double>& EMconstant = {},
                     const std::vector<double>& ECALgranularity = {},
                     const std::vector<double>& ECALEtaBoundaries = {},

                     const std::vector<double>& HADstochastic = {},
                     const std::vector<double>& HADconstant = {},
                     const std::vector<double>& HCALgranularity = {},
                     const std::vector<double>& HCALEtaBoundaries = {},

                     const std::vector<double>& CHlinear = {},
                     const std::vector<double>& CHconstant = {},
                     const std::vector<double>& CHMS = {},
                     const std::vector<double>& CHangular = {},
                     const std::vector<double>& trkEtaBoundaries = {},

                     unsigned maxReFit=50,
                     int verbose=0);

    arma::mat ptrans() const;
    arma::mat rawmat() const;
    double chisq() const;

    void killPU(arma::mat &ans); 
    void minimize();

private:
    void clear();
    void fillUncertainties();

    void doPrefit();

    bool clipValues();
    void initializeOptimizer();
    void buildLoss();
    
    jet recojet_, genjet_;

    arma::mat A_;
    std::vector<std::pair<unsigned, unsigned>> fitlocations_;

    double clipval_;

    std::shared_ptr<ParticleUncertainty> uncertainty_;
    std::shared_ptr<MatchingFilter> filter_;
    std::vector<std::shared_ptr<prefitter>> prefitters_;
    std::vector<bool> excludeGen_;

    bool recoverLostTracks_;

    unsigned maxReFit_;

    double PUexp_, PUpenalty_;
    std::shared_ptr<ChisqLossFCN> loss_;
    std::shared_ptr<MnMigrad> optimizer_;

    int verbose_;

    enum spatialLoss lossType_;
};

#endif
