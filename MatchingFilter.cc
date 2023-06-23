#include "SRothman/Matching/src/MatchingFilter.h"
#include "SRothman/SimonTools/src/deltaR.h"

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
