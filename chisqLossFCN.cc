#include "chisqLossFCN.h"
#include "matchingUtil.h"

double ChisqLossFCN::operator()(const std::vector<double>& data) const {
    //printf("top of loss\n");
    //fflush(stdout);
    arma::mat A = fullmat(recoPT.n_elem, genPT.n_elem, locations, data);
    //printf("1\n");
    //fflush(stdout);

    double lossPT = 0;
    double lossETA = 0;
    double lossPHI = 0;
    double lossPU = 0;

    arma::vec recoPT_pred = A * (genPT);
    //printf("2\n");
    //fflush(stdout);

    arma::vec PU = arma::conv_to<arma::vec>::from(recoPT_pred == 0);
    //printf("3\n");
    //fflush(stdout);

    arma::vec lossPTvec = (recoPT_pred - recoPT)/ errPT;
    lossPTvec %= (1-PU);
    lossPT = arma::dot(lossPTvec, lossPTvec);
    //printf("4\n");
    //fflush(stdout);

    recoPT_pred.replace(0, 1);
    //printf("5\n");
    //fflush(stdout);
    if (type == spatialLoss::TYPE1){
        arma::vec recoETA_pred = A * weightedGenETA;
        recoETA_pred /= recoPT_pred;
        arma::vec lossETAvec = (recoETA_pred - recoETA)/ errETA;
        lossETAvec %= (1-PU);
        lossETA = arma::dot(lossETAvec, lossETAvec);
        //printf("6\n");
        //fflush(stdout);

        arma::vec recoPHI_pred = A * weightedGenPHI;
        recoPHI_pred /= recoPT_pred;
        arma::vec lossPHIvec = (recoPHI_pred - recoPHI)/ errPHI;
        lossPHIvec %= (1-PU);
        lossPHI = arma::dot(lossPHIvec, lossPHIvec);
        //printf("7\n");
        //fflush(stdout);
    } else if (type == spatialLoss::TYPE2){
        arma::mat diffETA(NPReco, NPGen, arma::fill::none);
        diffETA.each_col() = recoETA;
        diffETA.each_row() -= arma::trans(genETA);
        lossETA = arma::accu(diffETA % diffETA);
        //printf("8\n");
        //fflush(stdout);

        arma::mat diffPHI(NPReco, NPGen, arma::fill::none);
        diffPHI.each_col() = recoPHI;
        diffPHI.each_row() -= arma::trans(genPHI);
        lossPHI = arma::accu(diffPHI % diffPHI);
        //printf("9\n");
        //fflush(stdout);
    }

    for(unsigned i=0; i<NPReco; ++i){
        //printf("PU loss for particle %u\n", i);
        //fflush(stdout);
        if(PU(i) == 0) continue;
        //printf("passed PU check\n");
        //fflush(stdout);
        //printf("ids[i] = %u\n", ids[i]);
        //fflush(stdout);
        double pt0 = PUpt0s[ids[i]];
        double exp = PUexps[ids[i]];
        double penalty = PUpenalties[ids[i]];
        //printf("got parameters\n");
        //fflush(stdout);

        lossPU += 2*exp*std::max(std::log(recoPT(i)/pt0), 0.0) + penalty;
    }
    //printf("10\n");
    //fflush(stdout);

    return lossPT + lossETA + lossPHI + lossPU;
};

