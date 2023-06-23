#include "matcher.h"

arma::mat matcher::ptrans(){
    arma::mat ans(A_);

    for(unsigned i=0; i<fitlocations_.size(); ++i){
        const auto& loc = fitlocations_[i];
        ans(loc.first, loc.second) = optimizer_->Value(i);
    }

    arma::vec colden = arma::sum(ans, 1);
    colden.replace(0, 1);
    ans.each_col() /= colden;

    arma::vec genpt = genjet_.ptvec();
    arma::vec recpt = recojet_.ptvec();
    genpt/=arma::accu(genpt);
    recpt/=arma::accu(recpt);

    arma::vec predpt = ans * genpt;

    for(unsigned iGen = 0; iGen < genjet_.nPart; ++iGen){
        for(unsigned iReco = 0; iReco < recojet_.nPart; ++iReco){
            if (predpt(iReco) > 0){
                ans(iReco, iGen) *= recpt(iReco)/predpt(iReco);
            }
        }
    }

    if(verbose_){
        printf("ptrans:\n");
        std::cout << ans;
        printf("GEN\n");
        std::cout << genjet_.ptvec().t();
        printf("ptrans * GEN\n");
        std::cout << (ans * genjet_.ptvec()).t();
        printf("RECO\n");
        std::cout << (recojet_.ptvec()).t();
    }
    return ans;
}

void matcher::killPU(arma::mat& ans){
    /*arma::rowvec denom = arma::sum(ans, 0);
    denom.replace(0, 1);
    arma::mat tmp(ans);
    tmp.each_row() /= denom;

    arma::vec recoPT_pred = tmp * realGenPT;
    arma::vec diff = arma::abs(recoPT_pred - realRecoPT)/realErrPT;
    arma::uvec bad = arma::find(diff > cutoff);
    if(verbose_){
        printf("killPU() found %llu bad reco particles\n", bad.n_elem);
        printf("absdif\n");
        std::cout << arma::trans(arma::abs(recoPT_pred - globalRecoPT)) << std::endl;
        printf("diff\n");
        std::cout << arma::trans(diff) << std::endl;
        printf("bad\n");
        std::cout << arma::trans(bad) << std::endl;
    }
    for(unsigned i=0; i<bad.n_elem; ++i){
        bool didAnything=false;
        unsigned iReco = bad(i);
        for(unsigned iGen = 0; iGen<ans.n_cols; ++iGen){
            if(ans(iReco, iGen)){
                ans(iReco, iGen) = 0;
                didAnything = true;
            }
        }
        if(verbose_ && didAnything){
            printf("killPU() killed reco %u\n", iReco);
        }
    }*/
}

const arma::mat matcher::A() const{
    arma::mat ans(A_);
    for(unsigned i=0; i<fitlocations_.size(); ++i){
        const auto& loc = fitlocations_[i];
        ans(loc.first, loc.second) = optimizer_->Value(i);
    }
    return ans;
}

void matcher::clear(){
    A_ = arma::mat(recojet_.particles.size(), 
                  genjet_.particles.size(), 
                  arma::fill::zeros);
}

void matcher::fillUncertainties(){
    for(particle& p : recojet_.particles){
        uncertainty_->addUncertainty(p, recojet_); 
    }
}

void matcher::doPrefit(const matcher* const previous){
                       
    std::vector<unsigned> floatingGen;
    if(previous){
        A_ = previous->A();
        for(unsigned iGen=0; iGen<genjet_.nPart; ++iGen){
            if(arma::accu(A_.col(iGen)) == 0){
                floatingGen.emplace_back(iGen);
            }
        }
    } else {
        floatingGen.resize(genjet_.nPart);
        std::iota(floatingGen.begin(), floatingGen.end(), 0);
    }

    fitlocations_.clear();
    for(unsigned iReco=0; iReco<recojet_.particles.size(); ++iReco){//foreach reco particle
        particle& reco = recojet_.particles[iReco];
        std::vector<unsigned> matchedgen;
        for(unsigned iGen : floatingGen){//foreach floating gen particle
            particle& gen = genjet_.particles[iGen];
            if(filter_->allowMatch(reco, gen, recojet_)){//if matching is allowed
                matchedgen.emplace_back(iGen);
            }//end if matching
        }//end foreach gen
        if(matchedgen.size() == 0){//if no gen particles match
            continue;
        } else if(matchedgen.size() == 1){//if exactly one match
            A_(iReco, matchedgen[0]) = 1;
        } else {//need to fit
            for(unsigned iGen : matchedgen){
                fitlocations_.emplace_back(iReco, iGen);
            }
        }//end switch(matchedgen.size())
    }//end foreach reco

    if(verbose_){
        if(previous){
            printf("from previous:\n");
            std::cout << previous->A() << std::endl;
        }
        printf("fixed by prefit:");
        std::cout << A_ << std::endl;
        arma::mat Q = arma::mat(A_.n_rows, A_.n_cols, arma::fill::zeros);
        for(auto& p : fitlocations_){
            Q(p.first, p.second) = 1;
        }
        printf("floating by prefit:");
        std::cout << Q << std::endl;
    }
}

void matcher::buildLoss(){
    if(fitlocations_.size()==0){
        return;
    }

    loss_ = std::make_unique<ChisqLossFCN>(
            A_, recojet_, genjet_, 
            fitlocations_, lossType_);
}

void matcher::initializeOptimizer(){
    if(!loss_){
        return;
    }
    MnUserParameters starting;
    char buffer[4];
    for(unsigned i=0; i<fitlocations_.size(); ++i){
        sprintf(buffer, "%u", i);
        starting.Add(buffer, 1.0, 1.0);
        starting.SetLowerLimit(buffer, 0.0);
    }
    optimizer_ = std::make_unique<MnMigrad>(*loss_, starting); 
}

void matcher::minimize(){
    if(!optimizer_){
        return;
    }
    unsigned iIter=0;
    do {
        (*optimizer_)();
    } while(clipValues() && ++iIter < maxReFit_-1); 
}

bool matcher::clipValues(){
    if(!optimizer_){
        return false;
    }
    bool didanything = false;
    arma::mat ans(A_);
    for(unsigned i=0; i<fitlocations_.size(); ++i){
        unsigned x=fitlocations_[i].first;
        unsigned y=fitlocations_[i].second;
        ans(x, y) = optimizer_->Value(i);
    }
    arma::vec colden = arma::sum(ans, 1);
    colden.replace(0, 1);
    ans.each_col() /= colden;

    for(unsigned i=0; i<fitlocations_.size(); ++i){
        unsigned x=fitlocations_[i].first;
        unsigned y=fitlocations_[i].second;
        double val = ans(x, y);
        if(val < clipval_){
            optimizer_->SetValue(i, 0);
            optimizer_->Fix(i);
            didanything = true;
        }
    }
    return didanything;
}
