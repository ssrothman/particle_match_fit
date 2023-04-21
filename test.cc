#include <stdio.h>
#include <memory>
#include <armadillo>
#include <random>
#include "chisqLossFCN.h"
#include "matcher.h"

int main(){
    printf("hello world\n");
    
    unsigned N = 8;

    auto j = std::make_shared<jet>();
    auto j_o = std::make_shared<jet>();

    gausJet(N, *j_o);
    auto ptrans = std::make_shared<arma::mat>(genJet(*j_o, *j, 
                    0.15, 0.05, 0.05,
                    0.20, 0.80, 0.00, 0.00, 0.2));

    std::uniform_real_distribution<double> unif(0.1, 2.0);
    std::default_random_engine re;

    std::vector<particle> recovec, genvec;
    for(unsigned i=0; i<j_o->nPart; ++i){
        recovec.emplace_back(j_o->pt[i], j_o->eta[i], j_o->phi[i],
                             0.15, 0.15, 0.15,
                             0, 0);
    }

    for(unsigned i=0; i<j->nPart; ++i){
        genvec.emplace_back(j->pt[i], j->eta[i], j->phi[i],
                            0, 0, 0,
                            0, 0);
    }

    printf("made vecs\n");
    matcher M(recovec, genvec);
    printf("true:\n");
    std::cout << *ptrans;
    printf("precomputed\n");
    std::cout << M.A;
    printf("to fit\n");
    std::cout << M.fitlocations;

}
