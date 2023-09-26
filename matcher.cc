#include "matcher.h"
#include "matchingUtil.h"

matcher::matcher(const jet& recojet,
                 const jet& genjet,

                 double clipval, 
   
                 const enum spatialLoss& loss,
                 const std::vector<double>& PUpt0s,
                 const std::vector<double>& PUexps,
                 const std::vector<double>& PUpenalties,

                 const std::string& uncertainty,

                 const std::vector<std::string>& filters,
                 const std::vector<double>& cutoffs, 
                 const std::vector<std::string>& prefitters,
                 const std::string& refiner,
                 const std::string& dropGenFilter,
                 const std::string& dropRecoFilter,

                 bool recoverLostTracks,
   
                 //uncertainty parameters
                 const std::vector<double>& EMstochastic, 
                 const std::vector<double>& EMnoise,
                 const std::vector<double>& EMconstant,
                 const std::vector<double>& ECALgranularityEta,
                 const std::vector<double>& ECALgranularityPhi,
                 const std::vector<double>& ECALEtaBoundaries,
   
                 const std::vector<double>& HADstochastic,
                 const std::vector<double>& HADconstant,
                 const std::vector<double>& HCALgranularityEta,
                 const std::vector<double>& HCALgranularityPhi,
                 const std::vector<double>& HCALEtaBoundaries,
   
                 const std::vector<double>& CHlinear,
                 const std::vector<double>& CHconstant,
                 const std::vector<double>& CHMSeta,
                 const std::vector<double>& CHMSphi,
                 const std::vector<double>& CHangularEta,
                 const std::vector<double>& CHangularPhi,
                 const std::vector<double>& trkEtaBoundaries,
   
                 unsigned maxReFit,
                 int verbose) :

            recojet_(recojet), genjet_(genjet),
            clipval_(clipval), 
            recoverLostTracks_(recoverLostTracks),
            maxReFit_(maxReFit), 
            PUpt0s_(PUpt0s),
            PUexps_(PUexps), PUpenalties_(PUpenalties),
            verbose_(verbose), lossType_(loss) {

    filters_ = std::make_unique<MatchingFilterEnsemble>(filters, cutoffs);

    uncertainty_ = ParticleUncertainty::get(
            uncertainty, 
            EMstochastic, EMnoise, EMconstant, 
            ECALgranularityEta, ECALgranularityPhi,
            ECALEtaBoundaries,

            HADstochastic, HADconstant, 
            HCALgranularityEta, HCALgranularityPhi,
            HCALEtaBoundaries,

            CHlinear, CHconstant, 
            CHMSeta, CHMSphi,
            CHangularEta, CHangularPhi,
            trkEtaBoundaries);

    prefitters_ = std::make_unique<PrefitterEnsemble>(prefitters, filters_);

    refiner_ = prefitRefiner::get(refiner);

    dropGenFilter_ = particleFilter::getParticleFilter(dropGenFilter);
    dropRecoFilter_ = particleFilter::getParticleFilter(dropRecoFilter);
    
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
        printf("rawmat * GEN\n");
        std::cout << (rawmat() * genjet_.ptvec()).t();
        printf("RECO\n");
        std::cout << (recojet_.ptvec()).t();
    }
    return ans;
}

void matcher::fillUncertainties(){
    if(verbose_>1){
        printf("\nUNCERTAINTIES:\n");
    }
    for(particle& p : recojet_.particles){
        uncertainty_->addUncertainty(p, recojet_); 
        if(verbose_>1){
            printf("(%0.5f, %0.5f, %0.5f)\n", p.dpt, p.deta, p.dphi);
        }
    }
    if(verbose_>1){
        printf("\n");
    }
}

void matcher::doPrefit(){
    std::vector<bool> usedGen(genjet_.particles.size(), false);

    std::unordered_map<unsigned, std::vector<unsigned>> recoToGen;

    for(unsigned iReco=0; iReco<recojet_.particles.size(); ++iReco){//foreach reco particle
        particle& reco = recojet_.particles[iReco];

        if(verbose_ > 9){
            printf("reco %u\n", iReco);
        }
        std::vector<unsigned> matchedgen = prefitters_->prefit(reco, genjet_);

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
            if(gen.charge==0 || gen.pdgid == 211){//missed muon tracks are gone for good
                continue;
            }

            if(verbose_ > 9){
                printf("recovering gen %u\n", iGen);
            }
            particle gencopy(gen);
            gencopy.charge = 0;

            for(unsigned iReco=0; iReco<recojet_.particles.size(); ++iReco){
                particle& reco = recojet_.particles[iReco];
                if(filters_->pass(reco, gencopy)){
                    if(verbose_>9){
                        printf("\treco %u: pass\n", iReco);
                    }
                    recoToGen[iReco].push_back(iGen);
                    usedGen[iGen] = true;
                } else {
                    if(verbose_>9){
                        printf("\treco %u: fail\n", iReco);
                    }
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
            PUpt0s_, PUexps_, PUpenalties_);
    if(verbose_>2){
        printf("loss mat\n");
        std::cout << rawmat();
        printf("fitlocations.size() = %lu\n", fitlocations_.size());
        fflush(stdout);
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
    greedyDropParticles(true);
    greedyDropParticles(false);

    iterativelyClip();
}

void matcher::greedyDropParticles(bool gen){
    unsigned maxI = gen ? genjet_.nPart : recojet_.nPart;

        
    for(unsigned i=0; i<maxI; ++i){
        if(gen && dropGenFilter_->pass(genjet_.particles[maxI-i-1])){
            testDrop(maxI-i-1, -1, false);
        } else if (!gen && dropRecoFilter_->pass(recojet_.particles[maxI-i-1])) {
            testDrop(-1, maxI-i-1, true);
        }
    }
}

void matcher::testDrop(int iGen, int iReco, bool allowInducedPU){
    if(!optimizer_){
        return;
    }

    if(iGen<0 && iReco<0){
        return;
    }

    arma::mat Araw_initial = rawmat();
    arma::vec PU_initial = arma::conv_to<arma::vec>::from(
            (Araw_initial * genjet_.ptvec()) == 0);

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

    arma::mat Araw_final = rawmat();
    arma::vec PU_final = arma::conv_to<arma::vec>::from(
            (Araw_final * genjet_.ptvec()) == 0);
    bool PUchanged = arma::any(PU_initial != PU_final);

    if(newchisq < savedchisq && (!PUchanged || allowInducedPU)){
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
