#include "MatchingFilter.h"
#include "matchingUtil.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/etaRegion.h"
#include "SRothman/SimonTools/src/util.h"

/*
 */

//#define DEBUG

class DRFilter : public MatchingFilter {
    public:
        DRFilter(const std::vector<double> etaBoundaries,
                 const std::vector<double> thresholds) : 
            MatchingFilter(),
            thresholds_(thresholds),
            etaBoundaries_(etaBoundaries) {
                /*printf("\nDRFilter: etaBoundaries_ = ");
                for (auto& eta : etaBoundaries_){
                    printf("%f, ", eta);
                }
                printf("\n");*/
            }
        ~DRFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            int region = getEtaRegion(reco.eta, etaBoundaries_);
            if(region < 0 || region > 2){
                printf("the bad region is %d\n", region);
                throw std::runtime_error("ThresholdFilter: region out of bounds");
            }
            double threshold = square(thresholds_[region]);
            return dR2(reco.eta, reco.phi, gen.eta, gen.phi) < threshold;
        }
    private:
        const std::vector<double> thresholds_;
        const std::vector<double> etaBoundaries_;
};

class ThresholdFilter : public MatchingFilter {
    public:
        ThresholdFilter(const std::vector<double>& etaBoundaries,
                        const std::vector<double>& thresholds) : 
            MatchingFilter(),
            thresholds_(thresholds),
            etaBoundaries_(etaBoundaries) {
                //printf("\nThresholdFilter: etaBoundaries_ = ");
                //for (auto& eta : etaBoundaries_){
                //    printf("%f, ", eta);
               // }
               // printf("\n");
            }
        ~ThresholdFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            int region = getEtaRegion(gen.eta, etaBoundaries_);
            if(region < 0 || region > 2){
                printf("the bad region is %d\n", region);
                throw std::runtime_error("ThresholdFilter: region out of bounds");
            }
            return reco.pt > thresholds_[region];
        }
    private:
        const std::vector<double> thresholds_;
        const std::vector<double> etaBoundaries_;
};

class ChargeSignFilter : public MatchingFilter {
    public:
        ChargeSignFilter() : MatchingFilter() {};
        ~ChargeSignFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.charge == reco.charge;
        }
};

class OppositeChargeSignFilter : public MatchingFilter {
    public:
        OppositeChargeSignFilter() : MatchingFilter() {};
        ~OppositeChargeSignFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.charge == -reco.charge;
        }
};

class ChargeMagnitudeFilter : public MatchingFilter {
    public:
        ChargeMagnitudeFilter() : MatchingFilter() {};
        ~ChargeMagnitudeFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return std::abs(gen.charge) == std::abs(reco.charge);
        }
};

class AnyChargedFilter : public MatchingFilter {
    public:
        AnyChargedFilter() : MatchingFilter() {};
        ~AnyChargedFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.charge != 0;
        }
};

class AnyChargedNoMuFilter : public MatchingFilter {
    public:
        AnyChargedNoMuFilter() : MatchingFilter() {};
        ~AnyChargedNoMuFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.charge != 0 && gen.pdgid != 13;
        }
};

class AnyNeutralFilter : public MatchingFilter {
    public:
        AnyNeutralFilter() : MatchingFilter() {};
        ~AnyNeutralFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.charge == 0;
        }
};

class AnyPhotonHardHadron1Filter : public MatchingFilter {
    public:
        AnyPhotonHardHadron1Filter() : MatchingFilter() {};
        ~AnyPhotonHardHadron1Filter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            if(gen.pdgid == 22){
                return true;
            } else if(gen.pdgid >= 100 && gen.charge==0 && gen.pt > 1){
                return true;
            } else {
                return false;
            }
        }
};

class AnyPhotonHardHadron2Filter : public MatchingFilter {
    public:
        AnyPhotonHardHadron2Filter() : MatchingFilter() {};
        ~AnyPhotonHardHadron2Filter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            if(gen.pdgid == 22){
                return true;
            } else if(gen.pdgid >= 100 && gen.charge==0 && gen.pt > 2){
                return true;
            } else {
                return false;
            }
        }
};

class AnyPhotonHardHadron5Filter : public MatchingFilter {
    public:
        AnyPhotonHardHadron5Filter() : MatchingFilter() {};
        ~AnyPhotonHardHadron5Filter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            if(gen.pdgid == 22){
                return true;
            } else if(gen.pdgid >= 100 && gen.charge==0 && gen.pt > 5){
                return true;
            } else {
                return false;
            }
        }
};

class AnyFilter : public MatchingFilter {
    public:
        AnyFilter() : MatchingFilter() {};
        ~AnyFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return true;
        }
};

class AnyPhotonFilter : public MatchingFilter {
    public:
        AnyPhotonFilter() : MatchingFilter() {};
        ~AnyPhotonFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.pdgid == 22;
        }
};

class AnyNeutralHadronFilter : public MatchingFilter {
    public:
        AnyNeutralHadronFilter() : MatchingFilter() {};
        ~AnyNeutralHadronFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.pdgid >= 100 && gen.charge == 0;
        }
};

class AnyChargedHadronFilter : public MatchingFilter {
    public:
        AnyChargedHadronFilter() : MatchingFilter() {};
        ~AnyChargedHadronFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.pdgid >= 100 && gen.charge != 0;
        }
};

class AnyHadronFilter : public MatchingFilter {
    public:
        AnyHadronFilter() : MatchingFilter() {};
        ~AnyHadronFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.pdgid >= 100;
        }
};

class AnyElectronFilter : public MatchingFilter {
    public:
        AnyElectronFilter() : MatchingFilter() {};
        ~AnyElectronFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.pdgid == 11;
        }
};

class AnyMuonFilter : public MatchingFilter {
    public:
        AnyMuonFilter() : MatchingFilter() {};
        ~AnyMuonFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return gen.pdgid == 13;
        }
};

std::shared_ptr<MatchingFilter> MatchingFilter::getFilter(const std::string& behavior){
    if (behavior == "ChargeSign"){
        return std::make_shared<ChargeSignFilter>();
    } else if (behavior == "OppositeChargeSign"){
        return std::make_shared<OppositeChargeSignFilter>();
    } else if (behavior == "ChargeMagnitude"){
        return std::make_shared<ChargeMagnitudeFilter>();
    } else if (behavior == "AnyCharged"){
        return std::make_shared<AnyChargedFilter>();
    } else if (behavior == "AnyChargedNoMu"){
        return std::make_shared<AnyChargedNoMuFilter>();
    } else if (behavior == "AnyNeutral"){
        return std::make_shared<AnyNeutralFilter>();
    } else if (behavior == "Any"){
        return std::make_shared<AnyFilter>();
    } else if (behavior == "AnyPhoton"){
        return std::make_shared<AnyPhotonFilter>();
    } else if (behavior == "AnyPhotonHardHadron1"){
        return std::make_shared<AnyPhotonHardHadron1Filter>();
    } else if (behavior == "AnyPhotonHardHadron2"){
        return std::make_shared<AnyPhotonHardHadron2Filter>();
    } else if (behavior == "AnyPhotonHardHadron5"){
        return std::make_shared<AnyPhotonHardHadron5Filter>();
    } else if (behavior == "AnyNeutralHadron"){
        return std::make_shared<AnyNeutralHadronFilter>();
    } else if (behavior == "AnyChargedHadron"){
        return std::make_shared<AnyChargedHadronFilter>();
    } else if (behavior == "AnyHadron"){
        return std::make_shared<AnyHadronFilter>();
    } else if (behavior == "AnyElectron"){
        return std::make_shared<AnyElectronFilter>();
    } else if (behavior == "AnyMuon"){
        return std::make_shared<AnyMuonFilter>();
    } else {
        printf("The bad behavior is %s\n", behavior.c_str());
        fflush(stdout);
        throw std::invalid_argument("Invalid matching filter type");
    }
}

MatchingFilterEnsemble::MatchingFilterEnsemble(
                const std::vector<std::string>& softflavorbehaviors,
                const std::vector<std::string>& hardflavorbehaviors,
                const std::vector<double>& behaviorthresholds,

                const std::vector<std::string>& chargebehaviors,

                const std::vector<double>& EM0thresholds,
                const std::vector<double>& HAD0thresholds,
                const std::vector<double>& HADCHthresholds,
                const std::vector<double>& ELEthresholds,
                const std::vector<double>& MUthresholds,

                const std::vector<double>& EM0dRcuts,
                const std::vector<double>& HAD0dRcuts,
                const std::vector<double>& HADCHdRcuts,
                const std::vector<double>& ELEdRcuts,
                const std::vector<double>& MUdRcuts,

                const std::vector<double>& ECALEtaBoundaries,
                const std::vector<double>& HCALEtaBoundaries,
                const std::vector<double>& trkEtaBoundaries) 
    : MatchingFilter(), 
      behaviorthresholds(behaviorthresholds){

    if(softflavorbehaviors.size() != 5 || 
       chargebehaviors.size() != 5 ||
       hardflavorbehaviors.size() != 5 ||
       behaviorthresholds.size() != 5){
        throw std::invalid_argument("Ensemble must have shape [EM0, HAD0, HADCH, ELE, MU]");
    }

    for(unsigned int i = 0; i < 5; i++){
        softflavorfilters.push_back(MatchingFilter::getFilter(
                    softflavorbehaviors[i]));
        hardflavorfilters.push_back(MatchingFilter::getFilter(
                    hardflavorbehaviors[i]));
        chargefilters.push_back(MatchingFilter::getFilter(chargebehaviors[i]));
    }
    dRfilters.push_back(std::make_shared<DRFilter>(ECALEtaBoundaries, EM0dRcuts));
    dRfilters.push_back(std::make_shared<DRFilter>(HCALEtaBoundaries, HAD0dRcuts));
    dRfilters.push_back(std::make_shared<DRFilter>(trkEtaBoundaries, HADCHdRcuts));
    dRfilters.push_back(std::make_shared<DRFilter>(trkEtaBoundaries, ELEdRcuts));
    dRfilters.push_back(std::make_shared<DRFilter>(trkEtaBoundaries, MUdRcuts));

    threshfilters.push_back(std::make_shared<ThresholdFilter>(
                ECALEtaBoundaries, EM0thresholds));
    threshfilters.push_back(std::make_shared<ThresholdFilter>(
                HCALEtaBoundaries, HAD0thresholds));
    threshfilters.push_back(std::make_shared<ThresholdFilter>(
                trkEtaBoundaries, HADCHthresholds));
    threshfilters.push_back(std::make_shared<ThresholdFilter>(
                trkEtaBoundaries, ELEthresholds));
    threshfilters.push_back(std::make_shared<ThresholdFilter>(
                trkEtaBoundaries, MUthresholds));
}

bool MatchingFilterEnsemble::pass(const particle& reco, const particle& gen){
    int index = pdgidToIndex(reco.pdgid, reco.charge);
    std::shared_ptr<MatchingFilter> flavorfilter;
    if(reco.pt < behaviorthresholds[index]){
        flavorfilter = softflavorfilters[index];
    } else{
        flavorfilter = hardflavorfilters[index];
    }
    return flavorfilter->pass(reco, gen) && 
           chargefilters[index]->pass(reco, gen) &&
           dRfilters[index]->pass(reco, gen) && 
           threshfilters[index]->pass(reco, gen);
}

int MatchingFilterEnsemble::pdgidToIndex(int pdgid, int charge) {
    if(pdgid==22){
        return 0;
    } else if (pdgid >=100 && charge==0){
        return 1;
    } else if (pdgid >=100 && charge!=0){
        return 2;
    } else if (pdgid == 11){
        return 3;
    } else if (pdgid == 13){
        return 4;
    } else {
        throw std::invalid_argument("Invalid particle type in MatchingFilterEnsemble");
    }
}
