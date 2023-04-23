#include <stdio.h>
#include <memory>
#include <armadillo>
#include <random>
#include "chisqLossFCN.h"
#include "simon_util_cpp/iterating.h"
#include "matcher.h"

int main(){
    printf("hello world\n");
    
    unsigned N = 10;

    auto j = std::make_shared<jet>();
    auto j_o = std::make_shared<jet>();

    gausJet(N, *j_o);
    auto ptrans = std::make_shared<arma::mat>(genJet(*j_o, *j, 
                    0.15, 0.05, 0.05,
                    0.10, 0.80, 0.10, 0.10, 0.2));

    std::uniform_real_distribution<double> unif(0.1, 2.0);
    std::default_random_engine re;

    std::vector<particle> recovec, genvec;
    arma::vec recoPT(j_o->nPart);
    for(unsigned i=0; i<j_o->nPart; ++i){
        recoPT(i) = j_o->particles[i].pt;
    }
    recovec = j_o->particles;

    arma::vec genPT(j->nPart);
    for(unsigned i=0; i<j->nPart; ++i){
        genPT(i) = j->particles[i].pt;
    }
    genvec = j->particles;


    printf("made vecs\n");
    matcher <spatialLoss::TYPE1> M1(recovec, genvec, 0.05, 2.0, true);
    matcher <spatialLoss::TYPE2> M2(recovec, genvec, 0.05, 2.0, true);
    printf("true:\n");
    std::cout << *ptrans;
    printf("precomputed\n");
    std::cout << M2.A;
    printf("recoToFit\n");
    printOrd(M1.recoToFit);
    printf("\ngenToFit\n");
    printOrd(M1.genToFit);
    printf("\n");

    //printf("preoptimization\n");
    //std::cout << M.optimizer->Parameters() << std::endl;
    printf("optimization..\n");
    M1.minimize();
    M2.minimize();
    //printf("postoptimization\n");
    //std::cout << M.optimizer->Parameters() << std::endl;

    printf("\nptrans1\n");
    std::cout << M1.ptrans();
    printf("\nptrans2\n");  
    std::cout << M2.ptrans();
    printf("genPT\n");
    std::cout << arma::trans(genPT);
    printf("recoPT\n");
    std::cout << arma::trans(recoPT);
    printf("A1 * gen\n");
    std::cout << arma::trans(M1.ptrans() * genPT);
    printf("A1 * 1\n");
    std::cout << arma::trans(M1.ptrans() * arma::ones(genPT.n_rows));

    arma::mat pt = M1.ptrans();
    printf("(gen, reco) = (%llu, %llu)\n", genPT.n_rows, recoPT.n_rows);
    printf("shape (%llu x %llu)\n", pt.n_rows, pt.n_cols);
    arma::rowvec rowsum = arma::sum(pt, 0);
    arma::vec colsum = arma::sum(pt, 1);
    printf("ROWSUM shape = (%llu, %llu)\n", rowsum.n_rows, rowsum.n_cols);
    printf("COLSUM shape = (%llu, %llu)\n", colsum.n_rows, colsum.n_cols);

    printf("sum(genPT): %0.3f\n", arma::sum(genPT));
    printf("sum(recoPT): %0.3f\n", arma::sum(recoPT));
    printf("sum(ptrans * genPT): %0.3f\n", arma::sum(pt * genPT));
}
