#include "SRothman/Matching/src/MatchingFilter.h"
#include "SRothman/SimonTools/src/deltaR.h"

bool MatchingFilter::passDR(const particle& part1, const particle& part2, const jet& j){
    double dR = std::sqrt(dR2(part1.eta,part1.phi, part2.eta,part2.phi));
    double err = std::sqrt(part1.deta*part1.deta + part1.dphi*part1.dphi);
    //printf("\tDR filter: dR: %0.3g, err: %0.3g, ratio: %0.3f\n", dR, err, dR/err);
    return dR/err < threshold_;
}

bool MatchingFilter::sameCharge(const particle& part1, const particle& part2, const jet& j){
    //printf("CHARGE filter: %i, %i\n", part1.charge, part2.charge);
    return std::abs(part1.charge) == std::abs(part2.charge);
}

bool MatchingFilter::sameChargeSign(const particle& part1, const particle& part2, const jet& j){
    //printf("CHARGESIGN filter: %i, %i\n", part1.charge, part2.charge);
    return part1.charge == part2.charge;
}

bool DRFilter::allowMatch(const particle& part1, const particle& part2, const jet& j){
    return passDR(part1, part2, j);
}

bool ChargeSignFilter::allowMatch(const particle& part1, const particle& part2, const jet& j){
    return passDR(part1, part2, j) && sameChargeSign(part1, part2, j);
}

bool ChargeFilter::allowMatch(const particle& part1, const particle& part2, const jet& j){
    return passDR(part1, part2, j) && sameCharge(part1, part2, j);
}

bool RealisticFilter::allowMatch(const particle& part1, const particle& part2, const jet& j){
    if (part1.pt < softPt_){
        return passDR(part1, part2, j);
    } else if (part1.pt < hardPt_){
        return passDR(part1, part2, j) && sameCharge(part1, part2, j);
    } else {
        return passDR(part1, part2, j);
    }
}
