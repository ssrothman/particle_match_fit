#include "matcher.h"

arma::mat matcher::ptrans(){
    arma::mat ans(A_);

    for(unsigned i=0; i<fitlocations_.n_rows; ++i){
        unsigned x=fitlocations_(i,0);
        unsigned y=fitlocations_(i,1);
        ans(x, y) = optimizer_->Value(i);
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
    for (unsigned i = 0; i < fitlocations_.n_rows; ++i){
        unsigned x = fitlocations_(i, 0);
        unsigned y = fitlocations_(i, 1);
        ans(x, y) = optimizer_->Value(i);
    }
    return ans;
}

void matcher::clear(){
    genToFit_.clear();
    recoToFit_.clear();
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
                       
    std::vector<unsigned>& floatingGen;
    std::vector<unsigned>& fixedGen;
    if(previous){
        A_ = previous->A();
        for(unsigned iGen=0; iGen<genjet_.nPart; ++iGen){
            if(arma::accu(A_.col(iGen)) == 0){
                floatingGen.emplace_back(iGen);
            } else {
                fixedGen.emplace_back(iGen);
            }
        }
    } else {
        floatingGen.resize(genjet_.nPart);
        std::iota(floatingGen.begin(), floatingGen.end(), 0);
    }

    std::unordered_set<unsigned> fittingReco, fittingGen;
    std::vector<std::pair<unsigned, unsigned>> toFit, fixedInFit;
    for(unsigned iReco=0; iReco<recojet_.particles.size(); ++iReco){//foreach reco particle
        particle& reco = recojet_.particles[iReco];
        std::vector<unsigned> matchedgen;
        for(unsigned iGen : floatingGen){//foreach floating gen particle
            particle& gen = genjet_.particles[iGen];
            if(filter_->pass(reco, gen)){//if matching is allowed
                matchedgen.emplace_back(iGen);
            }//end if matching
        }//end foreach gen
        if(matchedgen.size() == 0){//if no gen particles match
            continue;
        } else if(matchedgen.size() == 1){//if exactly one match
            A_(iReco, matchedgen[0]) = 1;
        } else {//need to fit
            for(unsigned iGen : matchedgen){
                toFit.emplace_back(iReco, iGen);
                fittingReco.insert(iReco);
                fittingGen.insert(iGen);
            }
            for(unsigned iGen : fixedGen){
                if(A_(iReco, iGen)){
                    fixedInFit.emplace_back(iReco, iGen);
                    fittingReco.insert(iReco);
                    fittingGen.insert(iGen);
                }
            }
        }//end switch(matchedgen.size())
    }//end foreach reco

    genToFit_.clear();
    genToFit_.reserve(fittingGen.size());
    genToFit_.insert(genToFit_.end(), fittingGen.begin(), 
                                      fittingGen.end());
    recoToFit_.clear();
    recoToFit_.reserve(fittingReco.size());
    recoToFit_.insert(recoToFit_.end(), fittingReco.begin(), 
                                        fittingReco.end());

    fitlocations_ = arma::umat(toFit.size(), 2, arma::fill::none);
    for(unsigned i=0; i<toFit.size(); ++i){
        fitlocations_(i, 0) = toFit[i].first;
        fitlocations_(i, 1) = toFit[i].second;
    }

    if(verbose_){
        if(previous){
            printf("from previous:\n");
            std::cout << previous->A() << std::endl;
        }
        printf("fixed by prefit:");
        std::cout << A_ << std::endl;
        arma::mat Q = arma::mat(A_.n_rows, A_.n_cols, arma::fill::zeros);
        for(auto& p : toFit){
            Q(p.first, p.second) = 1;
        }
        printf("floating by prefit:");
        std::cout << Q << std::endl;
    }
}

void matcher::buildLoss(){
    if(fitlocations_.n_rows==0){
        return;
    }
    std::unordered_map<unsigned, unsigned> recoIdxMap;
    arma::vec recoPT(recoToFit_.size(), arma::fill::none);
    arma::vec recoETA(recoToFit_.size(), arma::fill::none);
    arma::vec recoPHI(recoToFit_.size(), arma::fill::none);
    arma::vec errPT(recoToFit_.size(), arma::fill::none);
    arma::vec errETA(recoToFit_.size(), arma::fill::none);
    arma::vec errPHI(recoToFit_.size(), arma::fill::none);
    for(unsigned i=0; i<recoToFit_.size(); ++i){
        unsigned idx=recoToFit_[i];
        recoIdxMap[idx] = i;

        const particle& part = recojet_.particles[idx];
        recoPT[i] = part.pt;
        recoETA[i] = part.eta;
        recoPHI[i] = part.phi;
        errPT[i] = part.dpt;
        errETA[i] = part.deta;
        errPHI[i] = part.dphi;
    }

    std::unordered_map<unsigned, unsigned> genIdxMap;
    arma::vec genPT(genToFit_.size(), arma::fill::none);
    arma::vec genETA(genToFit_.size(), arma::fill::none);
    arma::vec genPHI(genToFit_.size(), arma::fill::none);
    for(unsigned i=0; i<genToFit_.size(); ++i){
        unsigned idx =genToFit_[i];
        genIdxMap[idx] = i;

        const particle& part = genjet_.particles[idx];
        genPT[i] = part.pt;
        genETA[i] = part.eta;
        genPHI[i] = part.phi;
    }

    arma::umat locations(fitlocations_.n_rows, 2u, arma::fill::none);
    for(unsigned i=0; i<fitlocations_.n_rows; ++i){
        locations(i, 0) = recoIdxMap[fitlocations_(i, 0)];
        locations(i, 1) = genIdxMap[fitlocations_(i, 1)];
    }

    loss_ = std::make_unique<ChisqLossFCN>(
        recoPT, recoETA, recoPHI,
        genPT, genETA, genPHI,
        errPT, errETA, errPHI,
        locations, lossType_);
}

void matcher::initializeOptimizer(){
    if(!loss_){
        return;
    }
    MnUserParameters starting;
    char buffer[4];
    for(unsigned i=0; i<fitlocations_.n_rows; ++i){
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
    for(unsigned i=0; i<fitlocations_.n_rows; ++i){
        unsigned x=fitlocations_(i,0);
        unsigned y=fitlocations_(i,1);
        ans(x, y) = optimizer_->Value(i);
    }
    arma::vec colden = arma::sum(ans, 1);
    colden.replace(0, 1);
    ans.each_col() /= colden;

    for(unsigned i=0; i<fitlocations_.n_rows; ++i){
        unsigned x=fitlocations_(i,0);
        unsigned y=fitlocations_(i,1);
        double val = ans(x, y);
        if(val < clipval_){
            optimizer_->SetValue(i, 0);
            optimizer_->Fix(i);
            didanything = true;
        }
    }
    return didanything;
}
