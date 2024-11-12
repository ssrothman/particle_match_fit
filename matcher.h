#ifndef PMF_MATCHER_H
#define PMF_MATCHER_H

#include <vector>
#include <unordered_set>
#include <iostream>

#include <string>

#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>

#include "SRothman/SimonTools/src/jets.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/util.h"

#include "chisqLossFCN.h"
#include "MatchingFilter.h"
#include "ParticleUncertainty.h"
#include "prefit.h"
#include "refinePrefit.h"
#include "particleFilter.h"

using namespace ROOT::Minuit2;

class matcher{
public:
    explicit matcher(const jet& recojet,
                     const jet& genjet,

                     double clipval, 

                     const enum spatialLoss& loss,
                     const std::vector<double>& PUpt0s,
                     const std::vector<double>& PUexps, 
                     const std::vector<double>& PUpenalties,

                     const std::string& uncertainty,

                     const std::vector<std::string>& softflavorfilters,
                     const std::vector<std::string>& hardflavorfilters,
                     const std::vector<double>& filterthresholds,

                     const std::vector<std::string>& chargefilters,

                     const std::vector<std::string>& dRfilters,

                     const std::vector<std::string>& prefitters,
                     const std::string& refiner,
                     const std::string& dropGenFilter,
                     const std::string& dropRecoFilter,

                     bool recoverLostTracks,
                     bool propagateLostTracks,
                     const std::vector<double>& HADCHrecoverThresholds,
                     const std::vector<double>& ELErecoverThresholds,
                     double Bz,

                     bool recoverLostHAD0,
                     const std::vector<double>& HAD0recoverThresholds,

                     const std::vector<double>& EMstochastic, 
                     const std::vector<double>& EMconstant,
                     const std::vector<double>& ECALgranularityEta,
                     const std::vector<double>& ECALgranularityPhi,
                     const std::vector<double>& ECALEtaBoundaries,

                     const std::vector<double>& HADstochastic,
                     const std::vector<double>& HADconstant,
                     const std::vector<double>& HCALgranularityEta,
                     const std::vector<double>& HCALgranularityPhi,
                     const std::vector<double>& HCALEtaBoundaries,

                     const std::vector<double>& CHlinear,
                     const std::vector<double>& CHconstant,
                     const std::vector<double>& CHMSeta,
                     const std::vector<double>& CHMSphi,
                     const std::vector<double>& CHangularEta,
                     const std::vector<double>& CHangularPhi,
                     const std::vector<double>& trkEtaBoundaries,

                     const std::vector<double>& EM0thresholds,
                     const std::vector<double>& HAD0thresholds,
                     const std::vector<double>& HADCHthresholds,
                     const std::vector<double>& ELEthresholds,
                     const std::vector<double>& MUthresholds,

                     const std::vector<double>& EM0constDR,
                     const std::vector<double>& EM0floatDR,
                     const std::vector<double>& EM0capDR,

                     const std::vector<double>& HAD0constDR,
                     const std::vector<double>& HAD0floatDR,
                     const std::vector<double>& HAD0capDR,

                     const std::vector<double>& HADCHconstDR,
                     const std::vector<double>& HADCHfloatDR,
                     const std::vector<double>& HADCHcapDR,

                     const std::vector<double>& ELEconstDR,
                     const std::vector<double>& ELEfloatDR,
                     const std::vector<double>& ELEcapDR,

                     const std::vector<double>& MUconstDR,
                     const std::vector<double>& MUfloatDR,
                     const std::vector<double>& MUcapDR,

                     unsigned maxReFit,
                     int verbose);

    Eigen::MatrixXd ptrans() const;
    Eigen::MatrixXd rawmat() const;
    double chisq() const;

    void minimize();

private:
    void fillUncertainties();

    void doPrefit();

    bool clipValues();
    void iterativelyClip();
    void initializeOptimizer();
    void buildLoss();

    void refineFit();
    void greedyDropParticles(bool gen);
    void testDrop(int iGen, int iReco, bool allowInducedPU);

    jet recojet_, genjet_;

    std::vector<std::pair<unsigned, unsigned>> fitlocations_;
    std::vector<bool> floating_;

    double clipval_;

    std::shared_ptr<ParticleUncertainty> uncertainty_;
    std::shared_ptr<MatchingFilterEnsemble> filters_;
    std::shared_ptr<PrefitterEnsemble> prefitters_;
    std::shared_ptr<prefitRefiner> refiner_;
    std::shared_ptr<particleFilter> dropGenFilter_;
    std::shared_ptr<particleFilter> dropRecoFilter_;

    std::vector<double> trkEtaBoundaries_;
    std::vector<double> ECALEtaBoundaries_;
    std::vector<double> HCALEtaBoundaries_;

    bool recoverLostTracks_;
    bool propagateLostTracks_;
    std::vector<double> HADCHrecoverThresholds_;
    std::vector<double> ELErecoverThresholds_;
    double Bz_;

    bool recoverLostHAD0_;
    std::vector<double> HAD0recoverThresholds_;

    unsigned maxReFit_;

    std::vector<double> PUpt0s_, PUexps_, PUpenalties_;
    std::shared_ptr<ChisqLossFCN> loss_;
    std::shared_ptr<MnMigrad> optimizer_;

    int verbose_;

    enum spatialLoss lossType_;
};

#endif
