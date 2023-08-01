#include "MatchingFilter.h"
#include "matchingUtil.h"
#include "SRothman/SimonTools/src/deltaR.h"

/*
 * Options: (enumerated in .cc rather than .h to avoid recompilation)
 *
 * DR: 
 */

static bool passDR(const particle& reco, const particle& gen, double threshold){
    return chisquared(reco, gen, false) < threshold;
}

static bool sameCharge(const particle& reco, const particle& gen){
    return reco.charge == gen.charge;
}

static bool sameChargeMagnitude(const particle& reco, const particle& gen){
    return std::abs(reco.charge) == std::abs(gen.charge);
}

class DRFilter : public MatchingFilter {
    public:
        DRFilter(double threshold) : MatchingFilter(threshold) {}
        ~DRFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return passDR(reco, gen, threshold_);
        }
};

class ChargeFilter : public MatchingFilter {
    public:
        ChargeFilter(double threshold) : MatchingFilter(threshold) {};
        ~ChargeFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return passDR(reco, gen, threshold_) && sameCharge(reco, gen);
        }
};

class ChargeMagnitudeFilter : public MatchingFilter {
    public:
        ChargeMagnitudeFilter(double threshold) : MatchingFilter(threshold) {};
        ~ChargeMagnitudeFilter() override {};
        bool pass(const particle& reco, const particle& gen) override{
            return passDR(reco, gen, threshold_) && sameChargeMagnitude(reco, gen);
        }
};

std::shared_ptr<MatchingFilter> MatchingFilter::getFilter(const std::string& behavior, double threshold){
    if(behavior == "DR"){
        return std::make_shared<DRFilter>(threshold);
    } else if (behavior == "Charge"){
        return std::make_shared<ChargeFilter>(threshold);
    } else if (behavior == "ChargeMagnitude"){
        return std::make_shared<ChargeMagnitudeFilter>(threshold);
    } else {
        throw std::invalid_argument("Invalid matching filter type");
    }
}

MatchingFilterEnsemble::MatchingFilterEnsemble(std::vector<std::string> behaviors, std::vector<double> thresholds){
    if(behaviors.size() != thresholds.size()){
        throw std::invalid_argument("Mismatched number of behaviors and thresholds");
    }

    if(behaviors.size() != 3){
        throw std::invalid_argument("Ensemble must have shape [EM0, HAD0, CH]");
    }

    for(unsigned int i = 0; i < behaviors.size(); i++){
        filters.push_back(MatchingFilter::getFilter(behaviors[i], thresholds[i]));
    }
}

bool MatchingFilterEnsemble::pass(const particle& reco, const particle& gen){
    if(reco.pdgid==22){
        return filters[0]->pass(reco, gen);
    } else if(reco.pdgid==130){
        return filters[1]->pass(reco, gen);
    } else if(reco.charge!=0) {
        return filters[2]->pass(reco, gen);
    } else {
        throw std::invalid_argument("Invalid particle type");
    }
}
