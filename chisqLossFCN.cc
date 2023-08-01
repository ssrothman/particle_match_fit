#include "chisqLossFCN.h"
#include "matchingUtil.h"

double ChisqLossFCN::operator()(const std::vector<double>& data) const {
    arma::mat A = fullmat(recoPT.n_elem, genPT.n_elem, locations, data);

    double lossPT = 0;
    double lossETA = 0;
    double lossPHI = 0;
    double lossPU = 0;

    arma::vec recoPT_pred = A * (genPT);

    arma::vec PU = arma::conv_to<arma::vec>::from(recoPT_pred == 0);

    arma::vec lossPTvec = (recoPT_pred - recoPT)/ errPT;
    lossPTvec %= (1-PU);
    lossPT = arma::dot(lossPTvec, lossPTvec);

    recoPT_pred.replace(0, 1);
    if (type == spatialLoss::TYPE1){
        arma::vec recoETA_pred = A * weightedGenETA;
        recoETA_pred /= recoPT_pred;
        arma::vec lossETAvec = (recoETA_pred - recoETA)/ errETA;
        lossETAvec %= (1-PU);
        lossETA = arma::dot(lossETAvec, lossETAvec);

        arma::vec recoPHI_pred = A * weightedGenPHI;
        recoPHI_pred /= recoPT_pred;
        arma::vec lossPHIvec = (recoPHI_pred - recoPHI)/ errPHI;
        lossPHIvec %= (1-PU);
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

    for(unsigned i=0; i<NPReco; ++i){
        if(PU(i) == 0) continue;
        double pt0 = PUpt0s[ids[i]];
        double exp = PUexps[ids[i]];
        double penalty = PUpenalties[ids[i]];

        lossPU += 2*exp*std::max(std::log(recoPT(i)/pt0), 0.0) + penalty;
    }

    return lossPT + lossETA + lossPHI + lossPU;
};

