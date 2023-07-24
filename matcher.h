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
#include "refinePrefit.h"

using namespace ROOT::Minuit2;

class matcher{
public:
    explicit matcher(const jet& recojet,
                     const jet& genjet,
                     const std::vector<bool>& excludeGen,
                     bool greedyDropMatches,
                     bool greedyDropGen,
                     bool greedyDropReco,

                     double clipval, 

                     const enum spatialLoss& loss,
                     const enum matchFilterType& filter,
                     const enum uncertaintyType& uncertainty,
                     const std::vector<enum prefitterType>& prefitters,
                     const enum prefitRefinerType& refiner,
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
    template<bool gen>
    void greedyDropParticles(){
        double bestchisq = chisq();

        unsigned nPart;
        if constexpr(gen){
            nPart = genjet_.nPart;
        } else {
            nPart = recojet_.nPart;
        }

        unsigned otherNPart;
        if constexpr(gen){
            otherNPart = recojet_.nPart;
        } else {
            otherNPart = genjet_.nPart;
        }

        for(unsigned iThis=0; iThis<nPart; ++iThis){//for each particle
            std::vector<double> oldstate = optimizer_->Params();
            arma::mat oldA(A_);

            bool anyFloating=false;
            std::vector<unsigned> changed;
            for(unsigned i=0; i<fitlocations_.size(); ++i){
                bool found;
                if constexpr(gen){
                    found = fitlocations_[i].second == iThis;
                } else {
                    found = fitlocations_[i].first == iThis;
                }
                if(found && optimizer_->Value(i) !=0){
                    anyFloating=true;
                    optimizer_->SetValue(i, 0);
                    optimizer_->Fix(i);
                    changed.emplace_back(i);
                }
            }

            for(unsigned iOther=0; iOther<otherNPart; ++iOther){
                double val;
                if constexpr(gen){
                    val = A_.at(iOther, iThis);
                } else {
                    val = A_.at(iThis, iOther);
                }
                if(val != 0){
                    if constexpr(gen){
                        A_.at(iOther, iThis) = 0;
                    } else {
                        A_.at(iThis, iOther) = 0;
                    }
                }
            }

            if(verbose_){
                printf("dropped particle %u:\n", iThis);
                std::cout << rawmat();
            }

            if(anyFloating){
                (*optimizer_)();
            }

            double newchisq = chisq();

            if(newchisq < bestchisq){//if improved
                if(verbose_){
                    printf("dropping particle %u reduced chisq from %f to %f\n", iThis, bestchisq, newchisq);
                }
                bestchisq = newchisq;
            } else {//else if got worse
                if(verbose_){
                    printf("dropping particle %u increased chisq from %f to %f\n", iThis, bestchisq, newchisq);
                }
                for(unsigned iMatch=0; iMatch < oldstate.size(); ++iMatch){
                    optimizer_->SetValue(iMatch, oldstate[iMatch]);
                }
                for(unsigned iMatch : changed){
                    optimizer_->Release(iMatch);
                }
                A_ = oldA;
            }//end if improved
        }//end for each particle
    }//end greedyDropParticles<gen>()


    void clear();
    void fillUncertainties();

    void doPrefit();

    bool clipValues();
    void initializeOptimizer();
    void buildLoss();

    void greedyDropMatches();
    
    jet recojet_, genjet_;

    arma::mat A_;
    std::vector<std::pair<unsigned, unsigned>> fitlocations_;

    double clipval_;

    std::shared_ptr<ParticleUncertainty> uncertainty_;
    std::shared_ptr<MatchingFilter> filter_;
    std::vector<std::shared_ptr<prefitter>> prefitters_;
    std::shared_ptr<prefitRefiner> refiner_;
    std::vector<bool> excludeGen_;
    bool greedyDropMatches_;
    bool greedyDropGen_;
    bool greedyDropReco_;

    bool recoverLostTracks_;

    unsigned maxReFit_;

    double PUexp_, PUpenalty_;
    std::shared_ptr<ChisqLossFCN> loss_;
    std::shared_ptr<MnMigrad> optimizer_;

    int verbose_;

    enum spatialLoss lossType_;
};

#endif
