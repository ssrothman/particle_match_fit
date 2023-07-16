#include "SRothman/Matching/src/MatchingFilter.h"
#include "SRothman/SimonTools/src/deltaR.h"

class DRFilter : public MatchingFilter {
    public:
        DRFilter(double threshold) : MatchingFilter(threshold) {}
        ~DRFilter() override {};
        bool allowMatch(const particle& reco, const particle& gen, const jet& j) override;
};

class ChargeSignFilter : public MatchingFilter {
    public:
        ChargeSignFilter(double threshold) : MatchingFilter(threshold) {};
        ~ChargeSignFilter() override {};
        bool allowMatch(const particle& reco, const particle& gen, const jet& j) override;
};

class ChargeFilter : public MatchingFilter {
    public:
        ChargeFilter(double threshold) : MatchingFilter(threshold) {};
        ~ChargeFilter() override {};
        bool allowMatch(const particle& reco, const particle& gen, const jet& j) override;
};

class RealisticFilter : public MatchingFilter {
    public:
        RealisticFilter(double threshold, double softPt, double hardPt) : 
            MatchingFilter(threshold),
            softPt_(softPt), hardPt_(hardPt) {};
        ~RealisticFilter() override {};
        bool allowMatch(const particle& reco, const particle& gen, const jet& j) override;
    private:
        double softPt_, hardPt_;
};

class LostTrackFilter : public MatchingFilter {
    public:
        LostTrackFilter(double threshold):
            MatchingFilter(threshold) {}
        ~LostTrackFilter() override {};
        bool allowMatch(const particle& reco, const particle& gen, const jet& j) override;
};

bool MatchingFilter::passDR(const particle& reco, const particle& gen, const jet& j){
    double dR = std::sqrt(dR2(reco.eta,reco.phi, gen.eta,gen.phi));
    double err = std::sqrt(reco.deta*reco.deta + reco.dphi*reco.dphi);
    //printf("\tDR filter: dR: %0.3g, err: %0.3g, ratio: %0.3f\n", dR, err, dR/err);
    return dR/err < threshold_;
}

bool MatchingFilter::sameCharge(const particle& reco, const particle& gen, const jet& j){
    //printf("CHARGE filter: %i, %i\n", reco.charge, gen.charge);
    return std::abs(reco.charge) == std::abs(gen.charge);
}

bool MatchingFilter::sameChargeSign(const particle& reco, const particle& gen, const jet& j){
    //printf("CHARGESIGN filter: %i, %i\n", reco.charge, gen.charge);
    return reco.charge == gen.charge;
}

bool DRFilter::allowMatch(const particle& reco, const particle& gen, const jet& j){
    return passDR(reco, gen, j);
}

bool ChargeSignFilter::allowMatch(const particle& reco, const particle& gen, const jet& j){
    return passDR(reco, gen, j) && sameChargeSign(reco, gen, j);
}

bool ChargeFilter::allowMatch(const particle& reco, const particle& gen, const jet& j){
    return passDR(reco, gen, j) && sameCharge(reco, gen, j);
}

bool RealisticFilter::allowMatch(const particle& reco, const particle& gen, const jet& j){
    if (reco.pt < softPt_){
        return passDR(reco, gen, j);
    } else if (reco.pt < hardPt_){
        return passDR(reco, gen, j) && sameCharge(reco, gen, j);
    } else {
        return passDR(reco, gen, j);
    }
}

bool LostTrackFilter::allowMatch(const particle& reco, const particle& gen, const jet& j){
    return passDR(reco, gen, j) && (sameCharge(reco, gen, j) || gen.charge==0);
}

std::shared_ptr<MatchingFilter> MatchingFilter::getFilter(const enum matchFilterType& behavior, double threshold){
    switch(behavior){
        case matchFilterType::DR:
            return std::make_shared<DRFilter>(threshold);
        case matchFilterType::CHARGESIGN:
            return std::make_shared<ChargeSignFilter>(threshold);
        case matchFilterType::CHARGE:
            return std::make_shared<ChargeFilter>(threshold);
        case matchFilterType::REALISTIC:
            throw std::invalid_argument("Realistic filter requires additional parameters");
        case matchFilterType::LOSTTRACK:
            return std::make_shared<LostTrackFilter>(threshold);
        default:
            throw std::invalid_argument("Invalid matching filter type");
    }
}

std::shared_ptr<MatchingFilter> MatchingFilter::getFilter(const enum matchFilterType& behavior, double threshold, double softPt, double hardPt){
    switch(behavior){
        case matchFilterType::DR:
            return std::make_shared<DRFilter>(threshold);
        case matchFilterType::CHARGESIGN:
            return std::make_shared<ChargeSignFilter>(threshold);
        case matchFilterType::CHARGE:
            return std::make_shared<ChargeFilter>(threshold);
        case matchFilterType::REALISTIC:
            return std::make_shared<RealisticFilter>(threshold, softPt, hardPt);
        case matchFilterType::LOSTTRACK:
            return std::make_shared<LostTrackFilter>(threshold);
        default:
            throw std::invalid_argument("Invalid matching filter type");
    }
}
