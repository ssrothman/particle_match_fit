#ifndef GENMATCHFCN_H
#define GENMATCHFN_H

#include "Minuit2/FCNBase.h"
#include <vector>
#include <memory>
#include <armadillo>

using vec_t = arma::vec;
using vecptr_t = std::shared_ptr<vec_t>;
using cvecptr_t = const std::shared_ptr<const vec_t>;

enum spatialLoss{
    TYPE1=0, //construct pT-weighted predicted pT, eta. Use in classic chisq loss
    TYPE2=1  //weight the chisq difference between recoETA, genETA by transfer matrix
};

template <enum spatialLoss type>
class ChisqLossFCN: public ROOT::Minuit2::FCNBase {
  private:
    //data defining fit problem
    cvecptr_t recoPT, recoETA, recoPHI;
    cvecptr_t genPT, genETA, genPHI;
    cvecptr_t errPT, errETA, errPHI;

    const size_t NPReco, NPGen;

    arma::vec weightedGenETA, weightedGenPHI;

    const double PUexp, PUpenalty;
    //constants for augmented lagrangian method
    //will be used to enforce sum_i Aij = 1
  public:
    explicit ChisqLossFCN(cvecptr_t recoPT, cvecptr_t recoETA, cvecptr_t recoPHI,
                         cvecptr_t genPT, cvecptr_t genETA, cvecptr_t genPHI,
                         cvecptr_t errPT, cvecptr_t errETA, cvecptr_t errPHI,
                         const double PUexp, const double PUpenalty):
      recoPT(recoPT), recoETA(recoETA), recoPHI(recoPHI),
      genPT(genPT), genETA(genETA), genPHI(genPHI),
      errPT(errPT), errETA(errETA), errPHI(errPHI),
      NPReco(recoPT->size()), NPGen(genPT->size()),
      weightedGenETA(*genPT % *genETA), 
      weightedGenPHI(*genPT % *genPHI),
      PUexp(PUexp), PUpenalty(PUpenalty) {}

    arma::mat vecToMat(const vec_t& data) const{
        arma::mat result(NPReco, NPGen, arma::fill::none);

        unsigned q=0;
        for(unsigned i=0; i<NPReco; ++i){
            for(unsigned j=0; j<NPGen; ++j){
                result(i, j) = data.at(q++);
            }
        }

        return result;
    } 

    double operator()(const std::vector<double>& data) const override{
        arma::mat A = vecToMat(data);
        arma::rowvec denom = arma::sum(A,0);
        denom.replace(0, 1);
        A.each_row() /= denom;

        double lossPT = 0;
        double lossETA = 0;
        double lossPHI = 0;
        double lossPU = 0;

        arma::vec recoPT_pred = A * (*genPT);
        arma::vec lossPTvec = (recoPT_pred - *recoPT)/ *errPT;
        lossPT = arma::dot(lossPTvec, lossPTvec);

        if constexpr(type == spatialLoss::TYPE1){
            arma::vec recoETA_pred = A * weightedGenETA;
            recoETA_pred /= recoPT_pred;
            arma::vec lossETAvec = (recoETA_pred - *recoETA)/ *errETA;
            lossETA = arma::dot(lossETAvec, lossETAvec);

            arma::vec recoPHI_pred = A * weightedGenPHI;
            recoPHI_pred /= recoPT_pred;
            arma::vec lossPHIvec = (recoPHI_pred - *recoPHI)/ *errPHI;
            lossPHI = arma::dot(lossPHIvec, lossPHIvec);
        } else if constexpr(type == spatialLoss::TYPE2){
            arma::mat diffETA(NPReco, NPGen, arma::fill::none);
            diffETA.each_col() = *recoETA;
            diffETA.each_row() -= arma::trans(*genETA);
            lossETA = arma::accu(diffETA % diffETA);

            arma::mat diffPHI(NPReco, NPGen, arma::fill::none);
            diffPHI.each_col() = *recoPHI;
            diffPHI.each_row() -= arma::trans(*genPHI);
            lossPHI = arma::accu(diffPHI % diffPHI);
        }

        lossPU = 0;

        return lossPT + lossETA + lossPHI + lossPU;
    };
    //error computation constant
    //should be 1.0 for our chisq likelihood 
    inline double Up() const override {return 1.0;}
};

#endif
