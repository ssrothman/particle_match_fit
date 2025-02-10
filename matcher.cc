#include "matcher.h"
#include "matchingUtil.h"
#include "SRothman/SimonTools/src/etaRegion.h"
#include "SRothman/SimonTools/src/isID.h"

#ifdef CMSSW
#include "CommonTools/BaseParticlePropagator/interface/RawParticle.h"
#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#endif

matcher::matcher(const simon::jet& recojet,
                 const simon::jet& genjet,

                 double clipval, 
   
                 const enum spatialLoss& loss,
                 const std::vector<double>& PUpt0s,
                 const std::vector<double>& PUexps,
                 const std::vector<double>& PUpenalties,

                 const std::string& uncertainty,

                 const std::vector<std::string>& softflavorfilters,
                 const std::vector<std::string>& hardflavorfilters,
                 const std::vector<double>& filterthresholds,

                 const std::vector<std::string>& chargefilters,

                 const std::vector<std::string>& dRfilters,

                 const std::vector<std::string>& prefitters,
                 const std::string& refiner,
                 const std::string& dropGenFilter,
                 const std::string& dropRecoFilter,

                 bool recoverLostTracks,
                 bool propagateLostTracks,
                 const std::vector<double>& HADCHrecoverThresholds,
                 const std::vector<double>& ELErecoverThresholds,
                 double Bz,

                 bool recoverLostHAD0,
                 const std::vector<double>& HAD0recoverThresholds,
   
                 //uncertainty parameters
                 const std::vector<double>& EMstochastic, 
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
   
                 const std::vector<double>& EM0thresholds,
                 const std::vector<double>& HAD0thresholds,
                 const std::vector<double>& HADCHthresholds,
                 const std::vector<double>& ELEthresholds,
                 const std::vector<double>& MUthresholds,

                 const std::vector<double>& EM0constDR,
                 const std::vector<double>& EM0floatDR,
                 const std::vector<double>& EM0capDR,

                 const std::vector<double>& HAD0constDR,
                 const std::vector<double>& HAD0floatDR,
                 const std::vector<double>& HAD0capDR,

                 const std::vector<double>& HADCHconstDR,
                 const std::vector<double>& HADCHfloatDR,
                 const std::vector<double>& HADCHcapDR,

                 const std::vector<double>& ELEconstDR,
                 const std::vector<double>& ELEfloatDR,
                 const std::vector<double>& ELEcapDR,

                 const std::vector<double>& MUconstDR,
                 const std::vector<double>& MUfloatDR,
                 const std::vector<double>& MUcapDR,

                 unsigned maxReFit,
                 int verbose) :

            recojet_(recojet), genjet_(genjet),
            clipval_(clipval), 
            trkEtaBoundaries_(trkEtaBoundaries),
            ECALEtaBoundaries_(ECALEtaBoundaries),
            HCALEtaBoundaries_(HCALEtaBoundaries),
            recoverLostTracks_(recoverLostTracks),
            propagateLostTracks_(propagateLostTracks),
            HADCHrecoverThresholds_(HADCHrecoverThresholds),
            ELErecoverThresholds_(ELErecoverThresholds),
            Bz_(Bz),
            recoverLostHAD0_(recoverLostHAD0),
            HAD0recoverThresholds_(HAD0recoverThresholds),
            maxReFit_(maxReFit), 
            PUpt0s_(PUpt0s),
            PUexps_(PUexps), PUpenalties_(PUpenalties),
            verbose_(verbose), lossType_(loss) {

    filters_ = std::make_unique<MatchingFilterEnsemble>(
            softflavorfilters,
            hardflavorfilters,
            filterthresholds,

            chargefilters,

            EM0thresholds, HAD0thresholds, 
            HADCHthresholds, ELEthresholds, MUthresholds,

            dRfilters,
            
            EM0constDR, EM0floatDR, EM0capDR,
            HAD0constDR, HAD0floatDR, HAD0capDR,
            HADCHconstDR, HADCHfloatDR, HADCHcapDR,
            ELEconstDR, ELEfloatDR, ELEcapDR,
            MUconstDR, MUfloatDR, MUcapDR,

            ECALEtaBoundaries, 
            HCALEtaBoundaries, 
            trkEtaBoundaries);

    uncertainty_ = ParticleUncertainty::get(
            uncertainty, 
            EMstochastic, EMconstant, 
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

Eigen::MatrixXd matcher::rawmat() const{
    if(optimizer_){
        return fullmat(recojet_.nPart, genjet_.nPart, 
                       fitlocations_, optimizer_->Params());
    } else {
        Eigen::MatrixXd result = fullmat(recojet_.nPart, genjet_.nPart, {}, {});
        return result;
    }
}

Eigen::MatrixXd matcher::ptrans() const {
    Eigen::MatrixXd ans = rawmat();

    Eigen::VectorXd genpt = genjet_.ptvec();
    Eigen::VectorXd recpt = recojet_.ptvec();
    genpt/=genjet_.sumpt;
    recpt/=recojet_.sumpt;

    Eigen::VectorXd predpt = ans * genpt;

    for(unsigned iGen = 0; iGen < genjet_.nPart; ++iGen){
        for(unsigned iReco = 0; iReco < recojet_.nPart; ++iReco){
            if (predpt(iReco) > 0){
                ans(iReco, iGen) *= recpt(iReco)/predpt(iReco);
            }
        }
    }

    if(verbose_ > 1){
        printf("ptrans:\n");
        std::cout << ans;
        printf("GEN\n");
        std::cout << genpt.transpose();
        printf("ptrans * GEN\n");
        std::cout << (ans * genpt).transpose();
        printf("rawmat * GEN\n");
        std::cout << (rawmat() * genpt).transpose();
        printf("RECO\n");
        std::cout << recpt.transpose();
    }
    return ans;
}

void matcher::fillUncertainties(){
    if(verbose_>1){
        printf("\nUNCERTAINTIES:\n");
    }
    for(simon::particle& p : recojet_.particles){
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
        simon::particle& reco = recojet_.particles[iReco];

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

    for(unsigned iGen=0; iGen<genjet_.particles.size(); ++iGen){
        if(usedGen[iGen]){
            continue;
        }
        const simon::particle& gen = genjet_.particles[iGen];

        if(isEM0(gen) || isMU(gen)){//can't recover photons or muons
            continue;
        }

        bool recover=false;
        if(recoverLostTracks_){
            int etaRegion = simon::getEtaRegion(gen.eta, trkEtaBoundaries_);
            if(etaRegion >= (int)trkEtaBoundaries_.size()-1){ //this can happen because genjets cluster in rapidity
                                                              //in this case the particle will be out of acceptance anyway
                                                              //so we can just skip it for recovery
                continue;
            }
            if(etaRegion < 0 || etaRegion >= (int)trkEtaBoundaries_.size()-1){
                printf("BIG PROBLEM: TRK ETA REGION OUT OF BOUNDS\n");
                printf("etaRegion: %d\n", etaRegion);
                printf("gen eta: %0.5f\n", gen.eta);
                printf("genjet eta: %0.5f\n", genjet_.eta);
                throw std::runtime_error("TRK ETA REGION OUT OF BOUNDS");
            }
            if(isELE(gen)){
                recover = recover || (gen.pt > ELErecoverThresholds_[etaRegion]);
            } else if(isHADCH(gen)){
                recover = recover || (gen.pt > HADCHrecoverThresholds_[etaRegion]);
            }
        }

        if(recoverLostHAD0_){
            if(isHAD0(gen)){
                int etaRegion = simon::getEtaRegion(gen.eta, HCALEtaBoundaries_);
                if(etaRegion < 0 || etaRegion >= (int)HCALEtaBoundaries_.size()-1){
                    printf("BIG PROBLEM: HCAL ETA REGION OUT OF BOUNDS\n");
                    printf("etaRegion: %d\n", etaRegion);
                    printf("gen eta: %0.5f\n", gen.eta);
                    printf("genjet eta: %0.5f\n", genjet_.eta);
                    throw std::runtime_error("HCAL ETA REGION OUT OF BOUNDS");
                }
                recover = recover || (gen.pt > HAD0recoverThresholds_[etaRegion]);
            }
        }

        if(!recover){
            continue;
        }

        if(verbose_ > 9){
            printf("recovering gen %u\n", iGen);
        }
        simon::particle gencopy(gen);
        if(isHAD0(gen)){
            gencopy.pdgid = 22;
        } else if(gencopy.charge != 0){

#ifdef CMSSW
            RawParticle tmppart;
            tmppart.setVertex(0, 0, 0, 0);
            double px = gen.pt * cos(gen.phi);
            double py = gen.pt * sin(gen.phi);
            double pz = gen.pt * sinh(gen.eta);
            double E = sqrt(px*px + py*py + pz*pz);
            tmppart.setMomentum(px, py, pz, E);
            tmppart.setCharge(gen.charge);

            BaseParticlePropagator prop(tmppart, 0, 0, 0);
            prop.setMagneticField(Bz_);

            gencopy.charge = 0;
            if(isELE(gen)){
                gencopy.pdgid = 22;
                if(propagateLostTracks_){
                    prop.propagateToEcalEntrance();
                }
            } else if(isHADCH(gen)){
                gencopy.pdgid = 130;
                if(propagateLostTracks_){
                    prop.propagateToHcalEntrance();
                }
            }
            if(propagateLostTracks_){
                gencopy.eta = prop.particle().eta();
                gencopy.phi = prop.particle().phi();
                gencopy.pt = prop.particle().pt();
            }
#endif
        }

        for(unsigned iReco=0; iReco<recojet_.particles.size(); ++iReco){
            simon::particle& reco = recojet_.particles[iReco];
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

    if(verbose_ > 1){
        Eigen::MatrixXd fixed(recojet_.nPart, genjet_.nPart);
        Eigen::MatrixXd floating(recojet_.nPart, genjet_.nPart);
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

    if(fitlocations_.size()==0){
        return;
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
    if(!loss_){
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
        printf("sum(gen) = %f\n", genjet_.ptvec().sum());
        printf("sum(A*gen) = %f\n", (rawmat()*genjet_.ptvec()).sum());
    }
}

bool matcher::clipValues(){
    if(!optimizer_){
        return false;
    }

    Eigen::MatrixXd A = rawmat();

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
    if(!loss_){
        return 0;
    }
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
    if(!optimizer_ || !loss_){
        return;
    }

    if(iGen<0 && iReco<0){
        return;
    }

    Eigen::MatrixXd Araw_initial = rawmat();
    Eigen::VectorXd PU_initial = ((Araw_initial * genjet_.ptvec()).array() == 0.).cast<double>();;

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

    Eigen::MatrixXd Araw_final = rawmat();
    Eigen::VectorXd PU_final = ((Araw_final * genjet_.ptvec()).array() == 0.0).cast<double>();
    bool PUchanged = (PU_initial.array() != PU_final.array()).any();

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
