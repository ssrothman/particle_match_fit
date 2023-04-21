#include <stdio.h>
#include "chisqLossFCN.h"
#include <memory>
#include <armadillo>

int main(){
    printf("hello world\n");
    
    unsigned NPReco = 10;
    unsigned NPGen = 13;
    double PUexp = 0;
    double PUpenalty = 0;

    auto recoPT = std::make_shared<arma::vec>(NPReco, arma::fill::randu);
    auto recoETA = std::make_shared<arma::vec>(NPReco, arma::fill::randu);
    auto recoPHI = std::make_shared<arma::vec>(NPReco, arma::fill::randu);

    auto errPT = std::make_shared<arma::vec>(NPReco, arma::fill::randu);
    auto errETA = std::make_shared<arma::vec>(NPReco, arma::fill::randu);
    auto errPHI = std::make_shared<arma::vec>(NPReco, arma::fill::randu);

    auto genPT = std::make_shared<arma::vec>(NPGen, arma::fill::randu);
    auto genETA = std::make_shared<arma::vec>(NPGen, arma::fill::randu);
    auto genPHI = std::make_shared<arma::vec>(NPGen, arma::fill::randu);

    std::vector<double> A(NPReco*NPGen, 1);

    ChisqLossFCN<spatialLoss::TYPE1> loss1(recoPT, recoETA, recoPHI,
                                           genPT, genETA, genPHI,
                                           errPT, errETA, errPHI,
                                           PUexp, PUpenalty);

    ChisqLossFCN<spatialLoss::TYPE2> loss2(recoPT, recoETA, recoPHI,
                                           genPT, genETA, genPHI,
                                           errPT, errETA, errPHI,
                                           PUexp, PUpenalty);

    printf("made loss\n");
    printf("loss1 = %0.3f\n", loss1(A));
    printf("loss2 = %0.3f\n", loss2(A));
}
