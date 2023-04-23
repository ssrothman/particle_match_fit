#ifndef PMF_MATCHER_H
#define PMF_MATCHER_H

#include <vector>
#include <unordered_set>
#include <iostream>

#include <armadillo>

#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>

#include "toyjets/common.h"

#include "simon_util_cpp/deltaR.h"
#include "simon_util_cpp/util.h"


using namespace ROOT::Minuit2;

template <enum spatialLoss type=TYPE1>
class matcher{
public:
    explicit matcher(const std::vector<particle>& recovec, 
                     const std::vector<particle>& genvec,
                     double clipval, double cutoff, bool matchCharge,
                     unsigned maxReFit=50) :
            A(recovec.size(), genvec.size(), arma::fill::zeros),
            globalGenPT(genvec.size()), globalRecoPT(recovec.size()),
            clipval(clipval), cutoff(cutoff), matchCharge(matchCharge),
            maxReFit(maxReFit){
        fitlocations = doPrefit(recovec, genvec);
        loss = buildLoss(recovec, genvec);
        optimizer = initializeOptimizer(recovec, genvec);

        for(unsigned i=0; i<genvec.size(); ++i){
            globalGenPT(i) = genvec[i].pt;
        }
        for(unsigned i=0; i<recovec.size(); ++i){
            globalRecoPT(i) = recovec[i].pt;
        }
    }

    arma::mat ptrans(){
        arma::mat ans(A);
        for(unsigned i=0; i<fitlocations.n_rows; ++i){
            unsigned x=fitlocations(i,0);
            unsigned y=fitlocations(i,1);
            ans(x, y) = optimizer->Value(i);
        }
        arma::vec colden = arma::sum(ans, 1);
        colden.replace(0, 1);
        ans.each_col() /= colden;
        ans.each_col() %= globalRecoPT;
        arma::rowvec rowden = arma::trans(globalGenPT)/arma::sum(globalGenPT);
        rowden.replace(0, 1);
        ans.each_row() /= rowden;
        return ans;
    }

    void minimize(){
        if(!optimizer){
            return;
        }
        unsigned iIter=0;
        do {
            (*optimizer)();
        } while(clipValues() && iIter++ < maxReFit); 
        printf("Took %u iterations\n", iIter);
    }

    bool clipValues(){
        if(!optimizer){
            return false;
        }

        bool didanything = false;
        for(unsigned i=0; i<fitlocations.n_rows; ++i){
            double val = optimizer->Value(i);
            if(val < clipval){
                optimizer->SetValue(i, 0);
                optimizer->Fix(i);
                didanything = true;
            }
        }
        return didanything;
    }

    std::unique_ptr<MnMigrad> initializeOptimizer(const std::vector<particle>& recovec,
                                                  const std::vector<particle>& genvec){
        if(!loss){
            return nullptr;
        }
        MnUserParameters starting;
        char buffer[4];
        for(unsigned i=0; i<fitlocations.n_rows; ++i){
            sprintf(buffer, "%u", i);
            starting.Add(buffer, 1.0, 1.0);
            starting.SetLowerLimit(buffer, 0.0);
        }
        return std::make_unique<MnMigrad>(*loss, starting); 
    }


    std::unique_ptr<ChisqLossFCN<type>> buildLoss(const std::vector<particle>& recovec,
                                 const std::vector<particle>& genvec){
        if(fitlocations.n_rows==0){
            return nullptr;
        }
        std::unordered_map<unsigned, unsigned> recoIdxMap;
        arma::vec recoPT(recoToFit.size(), arma::fill::none);
        arma::vec recoETA(recoToFit.size(), arma::fill::none);
        arma::vec recoPHI(recoToFit.size(), arma::fill::none);
        arma::vec errPT(recoToFit.size(), arma::fill::none);
        arma::vec errETA(recoToFit.size(), arma::fill::none);
        arma::vec errPHI(recoToFit.size(), arma::fill::none);
        for(unsigned i=0; i<recoToFit.size(); ++i){
            unsigned idx =recoToFit[i];
            recoIdxMap[idx] = i;

            const particle& part = recovec[idx];
            recoPT[i] = part.pt;
            recoETA[i] = part.eta;
            recoPHI[i] = part.phi;
            errPT[i] = part.dpt;
            errETA[i] = part.deta;
            errPHI[i] = part.dphi;

        }

        std::unordered_map<unsigned, unsigned> genIdxMap;
        arma::vec genPT(genToFit.size(), arma::fill::none);
        arma::vec genETA(genToFit.size(), arma::fill::none);
        arma::vec genPHI(genToFit.size(), arma::fill::none);
        for(unsigned i=0; i<genToFit.size(); ++i){
            unsigned idx =genToFit[i];
            genIdxMap[idx] = i;

            const particle& part = genvec[idx];
            genPT[i] = part.pt;
            genETA[i] = part.eta;
            genPHI[i] = part.phi;
        }

        arma::umat locations(fitlocations.n_rows, 2u, arma::fill::none);
        for(unsigned i=0; i<fitlocations.n_rows; ++i){
            locations(i, 0) = recoIdxMap[fitlocations(i, 0)];
            locations(i, 1) = genIdxMap[fitlocations(i, 1)];
        }

        return std::make_unique<ChisqLossFCN<type>>(
                                  recoPT, recoETA, recoPHI,
                                  genPT, genETA, genPHI,
                                  errPT, errETA, errPHI,
                                  0.0, 0.0,
                                  locations);
    }

    arma::umat doPrefit(const std::vector<particle>& recovec,
                        const std::vector<particle>& genvec){

        genToFit.clear();
        recoToFit.clear();
        A.fill(0);

        std::unordered_set<unsigned> genset;
        std::vector<std::pair<unsigned, unsigned>> locations;
        for(unsigned iReco=0; iReco<recovec.size(); ++iReco){
            std::vector<unsigned> matched = getMatched(recovec[iReco], genvec);
            if(matched.size()==0){
            } else if(matched.size()==1){
                A(iReco, matched[0]) = 1.0f;
            } else{
                //keep track of which gen&reco particles 
                //will need to be passed to the fitter
                recoToFit.emplace_back(iReco);
                genset.insert(matched.begin(), matched.end());

                //keep track of where the fit parameters 
                //will fit into the larger matrix
                for(unsigned i=0; i<matched.size(); ++i){
                    locations.emplace_back(iReco, matched[i]);
                }
            }
        }
        genToFit.insert(genToFit.end(), genset.begin(), genset.end());
        std::sort(genToFit.begin(), genToFit.end());

        arma::umat result(locations.size(), 2u, arma::fill::none);
        for(unsigned i=0; i<locations.size(); ++i){
            const auto& loc = locations[i];
            result(i, 0) = loc.first;
            result(i, 1) = loc.second;
        }
        return result;
    }

    std::vector<unsigned> getMatched(const particle& reco, 
                                     const std::vector<particle>& genvec) const {
        std::vector<unsigned> result;
        double dR2thresh = cutoff*(square(reco.deta)+square(reco.dphi));
        for(unsigned i=0; i<genvec.size(); ++i){
            const particle& gen = genvec[i];
            if(matchCharge && (gen.charge != reco.charge)){
                continue;
            }
            double dist2 = dR2(reco.eta, reco.phi, gen.eta, gen.phi);
            if(dist2 > dR2thresh){
                continue;
            }
            result.emplace_back(i);
        }
        return result;
    }

    arma::mat A;
    double clipval;
    double cutoff;
    bool matchCharge;
    
    unsigned maxReFit;

    arma::umat fitlocations;

    arma::vec globalGenPT, globalRecoPT;

    std::vector<unsigned> recoToFit;
    std::vector<unsigned> genToFit;

    std::unique_ptr<ChisqLossFCN<type>> loss;

    std::unique_ptr<MnMigrad> optimizer;
};

#endif
