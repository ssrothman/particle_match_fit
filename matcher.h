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

using namespace ROOT::Minuit2;

enum class matchFilterType{
    DR = 0,
    CHARGESIGN = 1,
    CHARGE = 2,
    REALISTIC = 3
};

enum class uncertaintyType{
    NAIVE = 0,
    STANDARD = 1,
    SMEAREDTRACKS = 2
};

class matcher{
public:
    explicit matcher(const jet& recojet,
                     const jet& genjet,

                     double clipval, 

                     enum spatialLoss loss,
                     enum matchFilterType filter,
                     enum uncertaintyType uncertainty,

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
                     int verbose=0,

                     const matcher *const previous=nullptr) :
            recojet_(recojet), genjet_(genjet),
            clipval_(clipval), maxReFit_(maxReFit), 
            verbose_(verbose), lossType_(loss) {

        if (filter == matchFilterType::DR){
            filter_ = std::make_unique<DRFilter>(cutoff);
        } else if(filter == matchFilterType::CHARGE){
            filter_ = std::make_unique<ChargeFilter>(cutoff);
        } else if(filter == matchFilterType::CHARGESIGN){
            filter_ = std::make_unique<ChargeSignFilter>(cutoff);
        } else if(filter == matchFilterType::REALISTIC){
            filter_ = std::make_unique<RealisticFilter>(cutoff, softPt, hardPt);
        } else {
            throw std::runtime_error("matcher: invalid filter type");
        }

        if(uncertainty == uncertaintyType::NAIVE){
            uncertainty_ = std::make_unique<NaiveParticleUncertainty>();
        } else if(uncertainty == uncertaintyType::STANDARD){
            uncertainty_ = std::make_unique<StandardParticleUncertainty>(
                EMstochastic,
                EMnoise,
                EMconstant,
                ECALgranularity,
                ECALEtaBoundaries,
                HADstochastic,
                HADconstant,
                HCALgranularity,
                HCALEtaBoundaries,
                CHlinear,
                CHconstant,
                CHMS,
                CHangular,
                trkEtaBoundaries);
        } else if (uncertainty == uncertaintyType::SMEAREDTRACKS){
            uncertainty_ = std::make_unique<StandardParticleUncertaintySmearedTracks>(
                EMstochastic,
                EMnoise,
                EMconstant,
                ECALgranularity,
                ECALEtaBoundaries,
                HADstochastic,
                HADconstant,
                HCALgranularity,
                HCALEtaBoundaries,
                CHlinear,
                CHconstant,
                CHMS,
                CHangular,
                trkEtaBoundaries);
        } else {
            throw std::runtime_error("matcher: invalid uncertainty type");
        }

        clear();
        fillUncertainties();
        doPrefit(previous);
        buildLoss();
        initializeOptimizer();
    }

    arma::mat ptrans();
    void killPU(arma::mat &ans); 
    void minimize();

    const arma::mat A() const;

private:
    void clear();
    void fillUncertainties();

    void doPrefit(const matcher *const previous=nullptr);

    bool clipValues();
    void initializeOptimizer();
    void buildLoss();
    void usePreviousFit(const matcher& previous);
    
    jet recojet_, genjet_;

    arma::mat A_;
    arma::umat fitlocations_;

    double clipval_;

    std::unique_ptr<ParticleUncertainty> uncertainty_;
    std::unique_ptr<MatchingFilter> filter_;

    unsigned maxReFit_;

    std::vector<unsigned> recoToFit_;
    std::vector<unsigned> genToFit_;

    std::unique_ptr<ChisqLossFCN> loss_;
    std::unique_ptr<MnMigrad> optimizer_;

    int verbose_;

    enum spatialLoss lossType_;


};

#endif
