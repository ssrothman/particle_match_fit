#include "MatchingFilter.h"
#include "matchingUtil.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/etaRegion.h"
#include "SRothman/SimonTools/src/util.h"
#include "SRothman/SimonTools/src/isID.h"

/*
 */

//#define DEBUG
#ifdef CMSSW_GIT_HASH

#endif

class FixedDRFilter : public MatchingFilter {
    public:
        FixedDRFilter(const std::vector<double>& etaBoundaries,
                      const std::vector<double>& thresholds) : 
            MatchingFilter(),
            thresholds_(thresholds),
            etaBoundaries_(etaBoundaries) {}
        ~FixedDRFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            int region = getEtaRegion(reco.eta, etaBoundaries_);
            if(region < 0 || region >= (int)(etaBoundaries_.size())-1){
                printf("the bad region is %d\n", region);
                throw std::runtime_error("FixedDRFilter: region out of bounds");
            }
            double threshold = square(thresholds_[region]);
            return dR2(reco.eta, reco.phi, gen.eta, gen.phi) < threshold;
        }
    private:
        const std::vector<double> thresholds_;
        const std::vector<double> etaBoundaries_;
};

class UncertaintyBasedDRFilter : public MatchingFilter {
    public:
        UncertaintyBasedDRFilter(const std::vector<double>& etaBoundaries,
                                 const std::vector<double>& thresholds) : 
            MatchingFilter(),
            thresholds_(thresholds),
            etaBoundaries_(etaBoundaries) {}
        ~UncertaintyBasedDRFilter() override {};

        bool pass(const particle& reco, const particle& gen) override{
            int region = getEtaRegion(reco.eta, etaBoundaries_);
            if(region < 0 || region >= (int)(etaBoundaries_.size()-1)){
                printf("the bad region is %d\n", region);
                throw std::runtime_error("FixedDRFilter: region out of bounds");
            }
            double threshold = square(thresholds_[region]);
            double unc = square(reco.deta) + square(reco.dphi);
            return dR2(reco.eta, reco.phi, gen.eta, gen.phi)/unc < threshold;
        }
    private:
        const std::vector<double> thresholds_;
        const std::vector<double> etaBoundaries_;
};

class TrackingDRFilter : public MatchingFilter {
    public:
        TrackingDRFilter(const std::vector<double>& etaBoundaries,
                         const std::vector<double>& hardTerm,
                         const std::vector<double>& invTerm,
                         const std::vector<double>& capTerm):
            MatchingFilter(),
            hardTerm_(hardTerm),
            invTerm_(invTerm),
            capTerm_(capTerm),
            etaBoundaries_(etaBoundaries) {}
        ~TrackingDRFilter() override {};

        bool pass(const particle& reco, const particle& gen) override {
            int region = getEtaRegion(reco.eta, etaBoundaries_);
            if(region < 0 || region >= (int)(etaBoundaries_.size())-1){
                printf("the bad region is %d\n", region);
                throw std::runtime_error("TrackingDRFilter: region out of bounds");
            }
            double hardTerm = hardTerm_[region];
            double invTerm = invTerm_[region];
            double capTerm = capTerm_[region];

            double threshold2 = square(hardTerm) + square(invTerm/gen.pt);
            threshold2 = std::min(threshold2, square(capTerm));

            return dR2(reco.eta, reco.phi, gen.eta, gen.phi) < threshold2;
        }
    private:
        const std::vector<double> hardTerm_;
        const std::vector<double> invTerm_;
        const std::vector<double> capTerm_;
        const std::vector<double> etaBoundaries_;
};

class ThresholdFilter : public MatchingFilter {
    public:
        ThresholdFilter(const std::vector<double>& etaBoundaries,
                        const std::vector<double>& thresholds) : 
            MatchingFilter(),
            thresholds_(thresholds),
            etaBoundaries_(etaBoundaries) {}
        ~ThresholdFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            int region = getEtaRegion(gen.eta, etaBoundaries_);
            if(region < 0 || region >= (int)(etaBoundaries_.size())-1){
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
            return gen.charge != 0 && !isMU(gen);
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
            if(isEM0(gen)){
                return true;
            } else if(isHAD0(gen) && gen.pt > 1){
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
            if(isEM0(gen)){
                return true;
            } else if(isHAD0(gen) && gen.pt > 2){
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
            if(isEM0(gen)){
                return true;
            } else if(isHAD0(gen) && gen.pt > 5){
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
            return isEM0(gen);
        }
};

class AnyNeutralHadronFilter : public MatchingFilter {
    public:
        AnyNeutralHadronFilter() : MatchingFilter() {};
        ~AnyNeutralHadronFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return isHAD0(gen);
        }
};

class AnyChargedHadronFilter : public MatchingFilter {
    public:
        AnyChargedHadronFilter() : MatchingFilter() {};
        ~AnyChargedHadronFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return isHADCH(gen);
        }
};

class AnyHadronFilter : public MatchingFilter {
    public:
        AnyHadronFilter() : MatchingFilter() {};
        ~AnyHadronFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return isHADCH(gen) || isHAD0(gen);
        }
};

class AnyElectronFilter : public MatchingFilter {
    public:
        AnyElectronFilter() : MatchingFilter() {};
        ~AnyElectronFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return isELE(gen);
        }
};

class AnyMuonFilter : public MatchingFilter {
    public:
        AnyMuonFilter() : MatchingFilter() {};
        ~AnyMuonFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return isMU(gen);
        }
};

std::shared_ptr<MatchingFilter> MatchingFilter::getDRFilter(
        const std::string& behavior,
        const std::vector<double>& etaBoundaries,
        const std::vector<double>& thresholds){
    if (behavior == "Fixed"){
        return std::make_shared<FixedDRFilter>(
                etaBoundaries, thresholds);
    } else if(behavior == "UncertaintyBased"){
        return std::make_shared<UncertaintyBasedDRFilter>(
                etaBoundaries, thresholds);
    } else if(behavior == "Tracking"){
        throw std::runtime_error("Tracking DR filter needs more arguments; use other getter");
    } else {
        throw std::runtime_error("Unknown DR filter");
    }
}

std::shared_ptr<MatchingFilter> MatchingFilter::getDRFilter(
        const std::string& behavior,
        const std::vector<double>& etaBoundaries,
        const std::vector<double>& hardTerm,
        const std::vector<double>& invTerm,
        const std::vector<double>& capTerm){
    if(behavior == "Tracking"){
        return std::make_shared<TrackingDRFilter>(
                etaBoundaries, hardTerm, invTerm, capTerm);
    } else {
        return getDRFilter(behavior, etaBoundaries, hardTerm);
    }
}

std::shared_ptr<MatchingFilter> MatchingFilter::getThresholdFilter(
        const std::vector<double>& boundaries,
        const std::vector<double>& thresholds){
    return std::make_shared<ThresholdFilter>(boundaries, thresholds);
}

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
                //flavor filters for soft particles
                const std::vector<std::string>& softflavorbehaviors,
                //flavor filters for hard particles
                const std::vector<std::string>& hardflavorbehaviors,
                //thresholds between "soft" and "hard" above
                const std::vector<double>& behaviorthresholds,

                //charge-matching behavior
                const std::vector<std::string>& chargebehaviors,

                //pT thresholds
                const std::vector<double>& EM0thresholds,
                const std::vector<double>& HAD0thresholds,
                const std::vector<double>& HADCHthresholds,
                const std::vector<double>& ELEthresholds,
                const std::vector<double>& MUthresholds,

                //dR cut behavior
                const std::vector<std::string>& DRbehaviors,
                //dR cut parameters
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

                //dR boundaries in ECAL, HCAL, tracker
                const std::vector<double>& ECALEtaBoundaries,
                const std::vector<double>& HCALEtaBoundaries,
                const std::vector<double>& trkEtaBoundaries) 
    : MatchingFilter(), 
      behaviorthresholds_(behaviorthresholds){

    if(softflavorbehaviors.size() != 5 || 
       chargebehaviors.size() != 5 ||
       hardflavorbehaviors.size() != 5 ||
       behaviorthresholds.size() != 5 ||
       DRbehaviors.size() != 5){
        throw std::invalid_argument("Ensemble must have shape [EM0, HAD0, HADCH, ELE, MU]");
    }

    for(unsigned int i = 0; i < 5; i++){
        softflavorfilters_.push_back(MatchingFilter::getFilter(
                    softflavorbehaviors[i]));
        hardflavorfilters_.push_back(MatchingFilter::getFilter(
                    hardflavorbehaviors[i]));
        chargefilters_.push_back(MatchingFilter::getFilter(chargebehaviors[i]));
    }

    std::vector<std::vector<double>> dRconst  = {
        EM0constDR, HAD0constDR, 
        HADCHconstDR, ELEconstDR, MUconstDR};
    std::vector<std::vector<double>> dRfloat  = {
        EM0floatDR, HAD0floatDR, 
        HADCHfloatDR, ELEfloatDR, MUfloatDR};
    std::vector<std::vector<double>> dRcap  = {
        EM0capDR, HAD0capDR, 
        HADCHcapDR, ELEcapDR, MUcapDR};

    std::vector<std::vector<double>> boundaries = {
        ECALEtaBoundaries, HCALEtaBoundaries, 
        trkEtaBoundaries, trkEtaBoundaries, trkEtaBoundaries};

    std::vector<std::vector<double>> thresholds = {
        EM0thresholds, HAD0thresholds, 
        HADCHthresholds, ELEthresholds, MUthresholds};

    for(unsigned int i = 0; i < 5; i++){
        dRfilters_.push_back(MatchingFilter::getDRFilter(
            DRbehaviors[i], boundaries[i],
            dRconst[i], dRfloat[i], dRcap[i]));

        threshfilters_.push_back(MatchingFilter::getThresholdFilter(boundaries[i], thresholds[i]));
    }
}

bool MatchingFilterEnsemble::pass(const particle& reco, const particle& gen){
    int index = pdgidToIndex(reco);
    std::shared_ptr<MatchingFilter> flavorfilter;
    if(reco.pt < behaviorthresholds_[index]){
        flavorfilter = softflavorfilters_[index];
    } else{
        flavorfilter = hardflavorfilters_[index];
    }
    return flavorfilter->pass(reco, gen) && 
           chargefilters_[index]->pass(reco, gen) &&
           dRfilters_[index]->pass(reco, gen) && 
           threshfilters_[index]->pass(reco, gen);
}

int MatchingFilterEnsemble::pdgidToIndex(const particle& p) {
    if(isEM0(p)){
        return 0;
    } else if (isHAD0(p)){
        return 1;
    } else if (isHADCH(p)){
        return 2;
    } else if (isELE(p)){
        return 3;
    } else if (isMU(p)){
        return 4;
    } else {
        throw std::invalid_argument("Invalid particle type in MatchingFilterEnsemble");
    }
}
