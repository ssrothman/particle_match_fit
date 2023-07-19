#include "matcher.h"

matcher::matcher(const jet& recojet,
                 const jet& genjet,
                 const std::vector<bool>& excludeGen,
   
                 double clipval, 
   
                 const  enum spatialLoss& loss,
                 const  enum matchFilterType& filter,
                 const  enum uncertaintyType& uncertainty,
                 const std::vector<enum prefitterType>& prefitters,
                 double PUexp, double PUpenalty,

                 bool recoverLostTracks,
   
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
                 int verbose) :

            recojet_(recojet), genjet_(genjet),
            clipval_(clipval), 
            excludeGen_(excludeGen),
            recoverLostTracks_(recoverLostTracks),
            maxReFit_(maxReFit), 
            PUexp_(PUexp), PUpenalty_(PUpenalty),
            verbose_(verbose), lossType_(loss) {

    filter_ = MatchingFilter::getFilter(filter, cutoff, softPt, hardPt);
    uncertainty_ = ParticleUncertainty::getUncertainty(uncertainty, 
                                                       EMstochastic, EMnoise, EMconstant, ECALgranularity, ECALEtaBoundaries,
                                                       HADstochastic, HADconstant, HCALgranularity, HCALEtaBoundaries,
                                                       CHlinear, CHconstant, CHMS, CHangular, trkEtaBoundaries);

    if(prefitters.size() !=3 ){
        throw std::runtime_error("matcher: invalid prefitter vector is wrong size");
    }

    prefitters_.resize(3);
    for(unsigned i=0; i<3; ++i){
        prefitters_[i] = prefitter::getPrefitter(prefitters[i], filter_, excludeGen);
    }
    
    clear();
    fillUncertainties();
    doPrefit();
    buildLoss();
    initializeOptimizer();
}

arma::mat matcher::rawmat() const{
    arma::mat ans(A_);

    for(unsigned i=0; i<fitlocations_.size(); ++i){
        const auto& loc = fitlocations_[i];
        ans(loc.first, loc.second) = optimizer_->Value(i);
    }

    arma::vec colden = arma::sum(ans, 1);
    colden.replace(0, 1);
    ans.each_col() /= colden;

    return ans;
}

arma::mat matcher::ptrans() const {
    arma::mat ans = rawmat();

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

void matcher::doPrefit(){
    fitlocations_.clear();

    std::vector<bool> usedGen(genjet_.particles.size(), false);

    for(unsigned iReco=0; iReco<recojet_.particles.size(); ++iReco){//foreach reco particle
        particle& reco = recojet_.particles[iReco];

        std::vector<unsigned> matchedgen;
        if(reco.pdgid==22){
            matchedgen = (*prefitters_[0])(reco, genjet_);
        }else if(reco.pdgid==130){
            matchedgen = (*prefitters_[1])(reco, genjet_);
        } else if(reco.charge!=0){
            matchedgen = (*prefitters_[2])(reco, genjet_);
        } else {
            throw std::runtime_error("matcher: invalid reco particle");
        }

        if(matchedgen.size() == 0){//if no gen particles match
            continue;
        } else if(matchedgen.size() == 1){//if exactly one match
            A_(iReco, matchedgen[0]) = 1;
            usedGen[matchedgen[0]] = true;
        } else {//need to fit
            for(unsigned iGen : matchedgen){
                fitlocations_.emplace_back(iReco, iGen);
                usedGen[iGen] = true;
            }
        }//end switch(matchedgen.size())
    }//end foreach reco

    if(recoverLostTracks_){
        for(unsigned iGen=0; iGen<genjet_.particles.size(); ++iGen){
            if(usedGen[iGen] || excludeGen_[iGen]){
                continue;
            }
            particle& gen = genjet_.particles[iGen];
            if(gen.charge==0){
                continue;
            }

            particle gencopy(gen);
            gencopy.charge = 0;

            for(unsigned iReco=0; iReco<recojet_.particles.size(); ++iReco){
                particle& reco = recojet_.particles[iReco];
                if(filter_->allowMatch(reco, gencopy, recojet_)){
                    fitlocations_.emplace_back(iReco, iGen);
                }
            }
        }
    }

    if(verbose_){
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

    loss_ = std::make_unique<ChisqLossFCN>(
            A_, recojet_, genjet_, 
            fitlocations_, lossType_,
            PUexp_, PUpenalty_);
    if(verbose_>2){
        printf("loss mat\n");
        std::cout << loss_->vecToMat(arma::vec(fitlocations_.size(), arma::fill::ones)) << std::endl;
        printf("loss = %f\n", loss_->operator()(std::vector<double>(fitlocations_.size(), 1.0)));
    }
}

void matcher::initializeOptimizer(){
    if(verbose_)
        printf("initializeOptimizer()\n");
    if(fitlocations_.size()==0){
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
    arma::mat ans = rawmat();

    for(unsigned i=0; i<fitlocations_.size(); ++i){
        unsigned x=fitlocations_[i].first;
        unsigned y=fitlocations_[i].second;
        double val = ans(x, y);
        if(val < clipval_ && val != 0){
            optimizer_->SetValue(i, 0);
            optimizer_->Fix(i);
            didanything = true;
        }
    }
    return didanything;
}

double matcher::chisq() const{
    if(!optimizer_){
        return loss_->operator()(std::vector<double>({}));
    }
    return (*loss_)(optimizer_->Params());
}
