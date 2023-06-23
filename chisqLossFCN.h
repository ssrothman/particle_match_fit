#ifndef PMF_CHISQ_LOSS_FCN_H
#define PMF_CHISQ_LOSS_FCN_H

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

    const arma::umat locations;

    const enum spatialLoss type;

  public:
    ChisqLossFCN() : 
        recoPT(1), recoETA(1), recoPHI(1),
        genPT(1), genETA(1), genPHI(1),
        errPT(1), errETA(1), errPHI(1),
        NPReco(1), NPGen(1),
        weightedGenETA(1), weightedGenPHI(1),
        locations(1, 1),
        type(TYPE1){}

    explicit ChisqLossFCN(const arma::vec& recoPT, 
                          const arma::vec& recoETA, 
                          const arma::vec& recoPHI,
                          const arma::vec& genPT, 
                          const arma::vec& genETA, 
                          const arma::vec& genPHI,
                          const arma::vec& errPT, 
                          const arma::vec& errETA,
                          const arma::vec& errPHI,
                          const arma::umat& locations,
                          const enum spatialLoss type):
      recoPT(recoPT), recoETA(recoETA), recoPHI(recoPHI),
      genPT(genPT), genETA(genETA), genPHI(genPHI),
      errPT(errPT), errETA(errETA), errPHI(errPHI),
      NPReco(recoPT.size()), NPGen(genPT.size()),
      weightedGenETA(genPT % genETA), 
      weightedGenPHI(genPT % genPHI),
      locations(locations),
      type(type){}

    arma::mat vecToMat(const arma::vec& data) const{
        arma::mat result(NPReco, NPGen, arma::fill::zeros);

        for(unsigned i=0; i<locations.n_rows; ++i){
            unsigned x=locations(i,0); 
            unsigned y=locations(i, 1);
            result(x, y) = data[i];
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

        arma::vec recoPT_pred = A * (genPT);
        arma::vec lossPTvec = (recoPT_pred - recoPT)/ errPT;
        lossPT = arma::dot(lossPTvec, lossPTvec);

        if (type == spatialLoss::TYPE1){
            arma::vec recoETA_pred = A * weightedGenETA;
            recoETA_pred /= recoPT_pred;
            arma::vec lossETAvec = (recoETA_pred - recoETA)/ errETA;
            lossETA = arma::dot(lossETAvec, lossETAvec);

            arma::vec recoPHI_pred = A * weightedGenPHI;
            recoPHI_pred /= recoPT_pred;
            arma::vec lossPHIvec = (recoPHI_pred - recoPHI)/ errPHI;
            lossPHI = arma::dot(lossPHIvec, lossPHIvec);
        } else if (type == spatialLoss::TYPE2){
            arma::mat diffETA(NPReco, NPGen, arma::fill::none);
            diffETA.each_col() = recoETA;
            diffETA.each_row() -= arma::trans(genETA);
            lossETA = arma::accu(diffETA % diffETA);

            arma::mat diffPHI(NPReco, NPGen, arma::fill::none);
            diffPHI.each_col() = recoPHI;
            diffPHI.each_row() -= arma::trans(genPHI);
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
