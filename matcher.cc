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

    arma::vec genpt = get_jet_pts(genjet_);
    arma::vec recpt = get_jet_pts(recojet_);
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
        std::cout << get_jet_pts(genjet_).t();
        printf("ptrans * GEN\n");
        std::cout << (ans * get_jet_pts(genjet_)).t();
        printf("RECO\n");
        std::cout << get_jet_pts(recojet_).t();
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
        0.0, 0.0,
        locations, lossType_);
}

void matcher::doPrefit(){
    genToFit_.clear();
    recoToFit_.clear();
    A_.fill(0);

    std::unordered_set<unsigned> genset;
    std::vector<std::pair<unsigned, unsigned>> locations;
    for(unsigned iReco=0; iReco<recojet_.particles.size(); ++iReco){
        std::vector<unsigned> matched = getMatched(recojet_.particles[iReco]);
        if(matched.size()==0){
        } else if(matched.size()==1){
            A_(iReco, matched[0]) = 1.0f;
        } else{
            //keep track of which gen&reco particles 
            //will need to be passed to the fitter
            recoToFit_.emplace_back(iReco);
            genset.insert(matched.begin(), matched.end());

            //keep track of where the fit parameters 
            //will fit into the larger matrix
            for(unsigned i=0; i<matched.size(); ++i){
                locations.emplace_back(iReco, matched[i]);
            }
        }
    }
    genToFit_.insert(genToFit_.end(), genset.begin(), genset.end());
    std::sort(genToFit_.begin(), genToFit_.end());

    fitlocations_=arma::umat(locations.size(), 2u, arma::fill::none);
    for(unsigned i=0; i<locations.size(); ++i){
        const auto& loc = locations[i];
        fitlocations_(i, 0) = loc.first;
        fitlocations_(i, 1) = loc.second;
    }
    if(verbose_){
        printf("\nprefit fixed matches:\n");
        std::cout << A_ << std::endl;
        printf("\nprefit floating matches:\n");
        arma::mat Q = arma::zeros<arma::mat>(recojet_.particles.size(), genjet_.particles.size());
        for(unsigned i=0; i<fitlocations_.n_rows; ++i){
            unsigned x=fitlocations_(i,0);
            unsigned y=fitlocations_(i,1);
            Q(x,y) = 1;
        }
        std::cout << Q << std::endl;
    }
}

std::vector<unsigned> matcher::getMatched(particle& reco){
    uncertainty_->addUncertainty(reco, recojet_);

    std::vector<unsigned> result;

    if (verbose_>1){
        printf("searching for matches for reco (%0.3f, %0.3f, %0.3f)\n", reco.pt, reco.eta, reco.phi);
    }
    for(unsigned i=0; i<genjet_.particles.size(); ++i){
        const particle& gen = genjet_.particles[i];
        if(filter_->allowMatch(reco, gen, recojet_)){
            result.emplace_back(i);
            if(verbose_>1){
                printf("\tPASS (%0.3f, %0.3f, %0.3f)\n", gen.pt, gen.eta, gen.phi);
            }
        } else {
            if(verbose_>1){
                printf("\tFAIL (%0.3f, %0.3f, %0.3f)\n", gen.pt, gen.eta, gen.phi);
            }
        }
    }
    return result;
}

