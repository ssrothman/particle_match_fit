#ifndef PMF_MATCHER_H
#define PMF_MATCHER_H

#include <vector>
#include <unordered_set>
#include <armadillo>
#include "toyjets/common.h"
#include "simon_util_cpp/deltaR.h"
#include "simon_util_cpp/util.h"

struct particle{
    double pt, eta, phi;
    double dpt, deta, dphi;
    unsigned pdgid; //absolute value
    int charge;

    particle(double pt, double eta, double phi,
             double dpt, double deta, double dphi,
             unsigned pdgid, int charge):
        pt(pt), eta(eta), phi(phi),
        dpt(dpt), deta(deta), dphi(dphi),
        pdgid(pdgid), charge(charge) {}
};

template <enum spatialLoss type=TYPE1>
class matcher{
public:
    explicit matcher(const std::vector<particle>& recovec, 
                     const std::vector<particle>& genvec) :
            A(recovec.size(), genvec.size(), arma::fill::zeros){

        fitlocations = doPrefit(recovec, genvec);
        loss = buildLoss(recovec, genvec);
    }

    void minimize();

    std::unique_ptr<ChisqLossFCN<type>> buildLoss(const std::vector<particle>& recovec,
                                 const std::vector<particle>& genvec){
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
        arma::vec genPT(recoToFit.size(), arma::fill::none);
        arma::vec genETA(recoToFit.size(), arma::fill::none);
        arma::vec genPHI(recoToFit.size(), arma::fill::none);
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
                recoToFit.reserve(matched.size());
                recoToFit.insert(recoToFit.end(), 1, 1);
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
        double dR2thresh = 2*(square(reco.deta)+square(reco.dphi));
        for(unsigned i=0; i<genvec.size(); ++i){
            const particle& gen = genvec[i];
            if(gen.charge != reco.charge){
                continue;
            }
            float dist2 = dR2(reco.eta, reco.phi, gen.eta, gen.phi);
            if(dist2 > dR2thresh){
                continue;
            }
            result.emplace_back(i);
        }
        return result;
    }

    arma::mat A;
    arma::umat fitlocations;

    std::vector<unsigned> recoToFit;
    std::vector<unsigned> genToFit;

    std::unique_ptr<ChisqLossFCN<type>> loss;
};

#endif
