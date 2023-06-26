#include "matcher.h"

matcher::matcher(const jet& recojet,
                 const jet& genjet,
   
                 double clipval, 
   
                 enum spatialLoss loss,
                 enum matchFilterType filter,
                 enum uncertaintyType uncertainty,
   
                 double cutoff, 
   
                 double softPt, double hardPt,
   
                 const std::vector<double>& EMstochastic, 
                 const std::vector<double>& EMnoise,
                 const std::vector<double>& EMconstant,
                 const std::vector<double>& ECALgranularity,
                 const std::vector<double>& ECALEtaBoundaries,
   
                 const std::vector<double>& HADstochastic,
                 const std::vector<double>& HADconstant,
                 const std::vector<double>& HCALgranularity,
                 const std::vector<double>& HCALEtaBoundaries,
   
                 const std::vector<double>& CHlinear,
                 const std::vector<double>& CHconstant,
                 const std::vector<double>& CHMS,
                 const std::vector<double>& CHangular,
                 const std::vector<double>& trkEtaBoundaries,
   
                 unsigned maxReFit,
                 int verbose,
   
                 const matcher *const previous) :
recojet_(recojet), genjet_(genjet),
            clipval_(clipval), maxReFit_(maxReFit), 
            verbose_(verbose), lossType_(loss) {

    if (filter == matchFilterType::DR){
        if(verbose_){
            printf("matcher: using DR filter with cutoff %f\n", cutoff);
        }
        filter_ = std::make_unique<DRFilter>(cutoff);
    } else if(filter == matchFilterType::CHARGE){
        if(verbose_){
            printf("matcher: using charge filter with cutoff %f\n", cutoff);
        }
        filter_ = std::make_unique<ChargeFilter>(cutoff);
    } else if(filter == matchFilterType::CHARGESIGN){
        if(verbose_){
            printf("matcher: using charge sign filter with cutoff %f\n", cutoff);
        }
        filter_ = std::make_unique<ChargeSignFilter>(cutoff);
    } else if(filter == matchFilterType::REALISTIC){
        if(verbose_){
            printf("matcher: using realistic filter with cutoff %f, softPt %f, hardPt %f\n", cutoff, softPt, hardPt);
        }
        filter_ = std::make_unique<RealisticFilter>(cutoff, softPt, hardPt);
    } else if(filter == matchFilterType::LOSTTRACK){
        if(verbose_){
            printf("matcher: using lost track filter with cutoff %f\n", cutoff);
        }
        filter_ = std::make_unique<LostTrackFilter>(cutoff);
    } else {
        throw std::runtime_error("matcher: invalid filter type");
    }

    if(uncertainty == uncertaintyType::NAIVE){
        if(verbose_){
            printf("matcher: using naive uncertainty\n");
        }
        uncertainty_ = std::make_unique<NaiveParticleUncertainty>();
    } else if(uncertainty == uncertaintyType::STANDARD){
        if(verbose_){
            printf("matcher: using standard uncertainty\n");
        }
        uncertainty_ = std::make_unique<StandardParticleUncertainty>(
            EMstochastic,
            EMnoise,
            EMconstant,
            ECALgranularity,
            ECALEtaBoundaries,
            HADstochastic,
            HADconstant,
            HCALgranularity,
            HCALEtaBoundaries,
            CHlinear,
            CHconstant,
            CHMS,
            CHangular,
            trkEtaBoundaries);
    } else if (uncertainty == uncertaintyType::SMEAREDTRACKS){
        if(verbose_){
            printf("matcher: using smeared tracks uncertainty\n");
        }
        uncertainty_ = std::make_unique<StandardParticleUncertaintySmearedTracks>(
            EMstochastic,
            EMnoise,
            EMconstant,
            ECALgranularity,
            ECALEtaBoundaries,
            HADstochastic,
            HADconstant,
            HCALgranularity,
            HCALEtaBoundaries,
            CHlinear,
            CHconstant,
            CHMS,
            CHangular,
            trkEtaBoundaries);
    } else {
        throw std::runtime_error("matcher: invalid uncertainty type");
    }

    clear();
    fillUncertainties();
    doPrefit(previous);
    buildLoss();
    initializeOptimizer();
}



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
    arma::rowvec denom = arma::sum(ans, 0);
    denom.replace(0, 1);
    ans.each_row() /= denom;
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
        printf("fixed by prefit:\n");
        std::cout << A_ << std::endl;
        arma::mat Q = arma::mat(A_.n_rows, A_.n_cols, arma::fill::zeros);
        for(auto& p : fitlocations_){
            Q(p.first, p.second) = 1;
        }
        printf("floating by prefit:\n");
        std::cout << Q << std::endl;
    }
}

void matcher::buildLoss(){
    if(verbose_){
        printf("buildLoss()\n");
    }
    if(fitlocations_.size()==0){
        return;
    }

    loss_ = std::make_unique<ChisqLossFCN>(
            A_, recojet_, genjet_, 
            fitlocations_, lossType_);
    if(verbose_>2){
        printf("loss mat\n");
        std::cout << loss_->vecToMat(arma::vec(fitlocations_.size(), arma::fill::ones)) << std::endl;
    }
}

void matcher::initializeOptimizer(){
    if(verbose_)
        printf("initializeOptimizer()\n");
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
    if(verbose_>2){
        printf("optimizer mat\n");
        std::cout << loss_->vecToMat(optimizer_->Params()) << std::endl;
    }
}

void matcher::minimize(){
    if(verbose_)
        printf("minimize()\n");
    if(!optimizer_){
        return;
    }
    unsigned iIter=0;
    do {
        (*optimizer_)();
        if(verbose_>2){
            printf("optimizer mat after iteration %u\n", iIter);
            std::cout << loss_->vecToMat(optimizer_->Params()) << std::endl;
        }
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
