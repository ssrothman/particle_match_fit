#include "matcher.h"
#include "matchingUtil.h"

matcher::matcher(const jet& recojet,
                 const jet& genjet,

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
        prefitters_[i] = prefitter::getPrefitter(prefitters[i], filter_);
    }

    refiner_ = prefitRefiner::getRefiner(refiner);
    
    fillUncertainties();
    doPrefit();
    buildLoss();
    initializeOptimizer();
}

arma::mat matcher::rawmat() const{
    if(optimizer_){
        return fullmat(recojet_.nPart, genjet_.nPart, 
                       fitlocations_, optimizer_->Params());
    } else {
        arma::mat result = fullmat(recojet_.nPart, genjet_.nPart, {}, {});
        return result;
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

void matcher::fillUncertainties(){
    for(particle& p : recojet_.particles){
        uncertainty_->addUncertainty(p, recojet_); 
    }
}

void matcher::doPrefit(){
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

        if(!matchedgen.empty()){
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
            const particle& gen = genjet_.particles[iGen];
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
        } else {
            for(unsigned iReco : p.second){
                fitlocations_.emplace_back(iReco, p.first);

                if(p.second.size()==1){
                    floating_.emplace_back(false);
                } else {
                    floating_.emplace_back(true);
                }
            }
        }
    }

    if(verbose_){
        arma::mat fixed(recojet_.nPart, genjet_.nPart, arma::fill::zeros);
        arma::mat floating(recojet_.nPart, genjet_.nPart, arma::fill::zeros);
        for(unsigned i=0; i<fitlocations_.size(); ++i){
            const auto& match = fitlocations_[i];
            if(floating_[i]){
                floating(match.first, match.second) = 1;
            } else {
                fixed(match.first, match.second) = 1;
            }
        }
        printf("fixed by prefit:\n");
        std::cout << fixed << std::endl;
        printf("floating by prefit:\n");
        std::cout << floating << std::endl;
    }
}

void matcher::buildLoss(){
    if(verbose_){
        printf("buildLoss()\n");
    }

    loss_ = std::make_unique<ChisqLossFCN>(
            recojet_, genjet_, 
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

    for(unsigned i=0; i<fitlocations_.size(); ++i){
        sprintf(buffer, "%u", i);
        if(floating_[i]){
            starting.Add(buffer, 0.5, 0.5, 0.0, 1.0);
            starting.Release(buffer);
        } else {
            starting.Add(buffer, 1.0);
            starting.Fix(buffer);
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

    (*optimizer_)();
    if(verbose_>1){
        printf("rough fit:\n");
        std::cout << rawmat();
    }
    
    refineFit();

    if (verbose_>1){
        printf("conservation of energy?\n");
        printf("sum(gen) = %f\n", arma::accu(genjet_.ptvec()));
        printf("sum(A*gen) = %f\n", arma::accu(rawmat()*genjet_.ptvec()));
    }
}

bool matcher::clipValues(){
    if(!optimizer_){
        return false;
    }

    arma::mat A = rawmat();

    bool didanything = false;
    for(unsigned i=0; i<fitlocations_.size(); ++i){
        const auto& match = fitlocations_[i];
        double val = A(match.first, match.second);
        if(floating_[i] && val < clipval_){
            optimizer_->SetValue(i, 0);
            optimizer_->Fix(i);
            floating_[i] = false;
            didanything = true;
        }
    }
    return didanything;
}

void matcher::iterativelyClip(){
    for(unsigned iIter=0; iIter < maxReFit_; ++iIter){
        if(!clipValues()){
            return;
        }
        if(verbose_>2){
            printf("optimizer mat after iteration %u\n", iIter);
            std::cout << rawmat();
        }
    }
}

double matcher::chisq() const{
    if(!optimizer_){
        return loss_->operator()(std::vector<double>({}));
    }
    return loss_->operator()(optimizer_->Params());
}

void matcher::refineFit(){
    if (greedyDropGen_){
        greedyDropParticles(true);
    } 

    if(greedyDropReco_){
        greedyDropParticles(false);
    }

    iterativelyClip();
}

void matcher::greedyDropParticles(bool gen){
    unsigned maxI = gen ? genjet_.nPart : recojet_.nPart;
    
    for(unsigned i=0; i<maxI; ++i){
        if(gen){
            testDrop(i, -1);
        } else {
            testDrop(-1, i);
        }
    }
}

void matcher::testDrop(int iGen, int iReco){
    if(!optimizer_){
        return;
    }

    if(iGen<0 && iReco<0){
        return;
    }

    MnUserParameters savedstate = optimizer_->Parameters();
    double savedchisq = chisq();
    std::vector<double> savedfloating(floating_.begin(), floating_.end());

    bool foundany=false;
    for(unsigned i=0; i<fitlocations_.size(); ++i){
        const auto& match = fitlocations_[i];
        if((int)match.second == iGen || (int)match.first == iReco){
            if(verbose_>2){
                printf("zeroing match %u -> %u\n", match.second, match.first);
            }
            optimizer_->RemoveLimits(i);
            optimizer_->SetValue(i, 0);
            optimizer_->Fix(i);
            floating_[i] = false;
            foundany = true;
        }
    }
    
    if(!foundany){
        return;
    }

    (*optimizer_)();
    
    double newchisq = chisq();

    if(newchisq < savedchisq){
        if(verbose_){
            printf("dropping gen %d, reco %d improved chisq from %f to %f\n", iGen, iReco, savedchisq, newchisq);
            std::cout << rawmat();
        }
    } else {
        if(verbose_){
            printf("dropping gen %d, reco %d worsened chisq from %f to %f\n", iGen, iReco, savedchisq, newchisq);
            std::cout << rawmat();
        }
        optimizer_ = std::make_unique<MnMigrad>(*loss_, savedstate);
        for(unsigned i=0; i<floating_.size(); ++i){
            floating_[i] = savedfloating[i];
        }
    }
}
