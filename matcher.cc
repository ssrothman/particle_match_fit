#include "matcher.h"
#include "matchingUtil.h"

/*
 * NOTES FOR TOMORROW:
 *
 * greedy things are fucked up
 *      fixing/unfixing parameters is sketchy
 *      changes to base A don't propagate to minimization problem
 * 
 * prefit refinement is fucked up
 *      just reutns a bunch of nothing
 *      no idea why at the moment
 */

matcher::matcher(const jet& recojet,
                 const jet& genjet,
                 const std::vector<bool>& excludeGen,
                 bool greedyDropMatches,
                 bool greedyDropGen,
                 bool greedyDropReco,
   
                 double clipval, 
   
                 const enum spatialLoss& loss,
                 const enum matchFilterType& filter,
                 const enum uncertaintyType& uncertainty,
                 const std::vector<enum prefitterType>& prefitters,
                 const enum prefitRefinerType& refiner,

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
            greedyDropMatches_(greedyDropMatches),
            greedyDropGen_(greedyDropGen),
            greedyDropReco_(greedyDropReco),
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

    refiner_ = prefitRefiner::getRefiner(refiner);
    
    clear();
    fillUncertainties();
    doPrefit();
    buildLoss();
    initializeOptimizer();
}

arma::mat matcher::rawmat() const{
    if(optimizer_){
        return fullmat(A_, fitlocations_, optimizer_->Params());
    } else {
        return fullmat(A_, {}, {});
    }
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

    std::unordered_map<unsigned, std::vector<unsigned>> recoToGen;

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

        if(matchedgen.size()){
            recoToGen[iReco] = matchedgen;
            for(unsigned iGen : matchedgen){
                usedGen[iGen] = true;
            }
        }
    }//end foreach reco

    if(recoverLostTracks_){
        for(unsigned iGen=0; iGen<genjet_.particles.size(); ++iGen){
            if(usedGen[iGen]){
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
                    recoToGen[iReco].push_back(iGen);
                    usedGen[iGen] = true;
                }
            }
        }
    }

    prefitRefiner::matchMap genToReco = refiner_->refine(recoToGen, recojet_, genjet_);

    for(const auto& p : genToReco){
        if(p.second.size()==0){
            continue;
        //} else if(p.second.size()==1){
        //    //A_(p.second[0], p.first) = 1;
        } else {
            for(unsigned iReco : p.second){
                fitlocations_.emplace_back(iReco, p.first);
            }
        }
    }

    if(verbose_){
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

    loss_ = std::make_unique<ChisqLossFCN>(
            A_, recojet_, genjet_, 
            fitlocations_, lossType_,
            PUexp_, PUpenalty_);
    if(verbose_>2){
        printf("loss mat\n");
        std::cout << rawmat();
        printf("loss = %f\n", loss_->operator()(std::vector<double>(fitlocations_.size(), 1.0)));
    }
}

void matcher::initializeOptimizer(){
    if(verbose_>1)
        printf("initializeOptimizer()\n");
    if(fitlocations_.size()==0){
        return;
    }
    MnUserParameters starting;
    char buffer[4];
    std::unordered_map<unsigned, std::vector<unsigned>> genToMatchIdx;
        ;
    for(unsigned i=0; i<fitlocations_.size(); ++i){
        sprintf(buffer, "%u", i);
        starting.Add(buffer, 0.5, 0.5, 0.0, 1.0);
        genToMatchIdx[fitlocations_[i].second].emplace_back(i);
    }
    for(const auto& match : genToMatchIdx){
        if(match.second.size()==1){
            unsigned iMatch = match.second[0];
            starting.RemoveLimits(iMatch);
            starting.SetValue(iMatch, 1.0);
            starting.Fix(iMatch);
        }
    }
    optimizer_ = std::make_unique<MnMigrad>(*loss_, starting); 
    if(verbose_>2){
        printf("optimizer mat\n");
        std::cout << rawmat();
    }
}

void matcher::minimize(){
    if(verbose_>1)
        printf("minimize()\n");
    if(!optimizer_){
        return;
    }
    unsigned iIter=0;
    do {
        (*optimizer_)();
        if(verbose_>2){
            printf("optimizer mat after iteration %u\n", iIter);
            std::cout << rawmat();
        }
    } while(clipValues() && ++iIter < maxReFit_-1); 

    if(greedyDropGen_){
        greedyDropParticles<true>();
    }

    if(greedyDropMatches_){
        greedyDropMatches();
    }

    if(greedyDropReco_){
        greedyDropParticles<false>();
    }

    if (verbose_>1){
        printf("conservation of energy?\n");
        printf("sum(gen) = %f\n", arma::accu(genjet_.ptvec()));
        printf("sum(A*gen) = %f\n", arma::accu(rawmat()*genjet_.ptvec()));
    }
}

void matcher::greedyDropMatches(){
    double bestchisq = chisq();
    for(unsigned iMatch=0; iMatch < fitlocations_.size(); ++iMatch){//for each match
        double savedval = optimizer_->Value(iMatch);
        if(savedval == 0){
            continue;
        }
        optimizer_->SetValue(iMatch, 0);
        optimizer_->Fix(iMatch);
        (*optimizer_)();
        double newchisq = chisq();
        if(verbose_>2){

            printf("dropping match (%u, %u):\n", fitlocations_[iMatch].first, fitlocations_[iMatch].second);
            std::cout << rawmat() << std::endl;
        }
        if(newchisq < bestchisq){
            if(verbose_>2){
                printf("\treduced chisq from %f to %f\n", bestchisq, newchisq);
            }
            bestchisq = newchisq;
        } else {
            if(verbose_>2){
                printf("\tincreased chisq from %f to %f\n", bestchisq, newchisq);
            }
            optimizer_->SetValue(iMatch, savedval);
            optimizer_->Release(iMatch);
        }
    }
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
    return loss_->operator()(optimizer_->Params());
}
