#include "ParticleUncertainty.h"
#include "SRothman/SimonTools/src/util.h"

static int getEtaRegion(const double& eta, const std::vector<double>& boundaries){
    for(unsigned i=0; i<boundaries.size(); ++i){
        if(std::abs(eta) < boundaries[i]){
            return i-1;
        }
    }
    return boundaries.size();
}

static double caloResolution(const double& eta, const double& pt, const double& A, const double& B, const double& C){
    double E = pt * std::cosh(eta);
    return std::sqrt(square(A/std::sqrt(E)) + square(B/E) + square(C)) * pt;
}

static double trkResolution(const double& pt, const double& A, const double& B){
    return std::sqrt(square(A*pt) + square(B)) * pt;
}

static double trkAngularResolution(const double& pt, const double& A, const double& B){
    return std::sqrt(square(A/pt) + square(B));
}

void NaiveParticleUncertainty::addUncertainty(particle& part, const jet& j){
    part.dpt = 0.1 * part.pt;
    part.dphi = 0.05;
    part.deta = 0.05;
}

void RealisticParticleUncertainty::addUncertaintyToNeutralEM(particle& part){
    int region = getEtaRegion(part.eta, ECALEtaBoundaries_);

    part.dpt = caloResolution(part.eta, part.pt, EMstochastic_[region], EMnoise_[region], EMconstant_[region]);
    part.dphi = ECALgranularity_[region];
    part.deta = part.dphi;
}

void RealisticParticleUncertainty::addUncertaintyToNeutralHadron(particle& part){
    int region = getEtaRegion(part.eta, HCALEtaBoundaries_);

    double E = part.pt * std::cosh(part.eta);
    part.dpt = caloResolution(part.eta, part.pt, HADstochastic_[region], 0, HADconstant_[region]);
    part.dphi = HCALgranularity_[region];
    part.deta = part.dphi;
}

void RealisticParticleUncertainty::addUncertaintyToCharged(particle & part){
    int region = getEtaRegion(part.eta, trkEtaBoundaries_);

    double trkRes = trkResolution(part.pt, CHlinear_[region], CHconstant_[region]);
    double caloRes = 9e9;
    if(part.pdgid==11){
        caloRes = caloResolution(part.eta, part.pt, 
                                        EMstochastic_[region], 
                                        EMnoise_[region], 
                                        EMconstant_[region]);
    } else if(part.pdgid!=13){//muons are tracker-only objects
        caloRes = caloResolution(part.eta, part.pt, 
                                HADstochastic_[region], 
                                0, 
                                HADconstant_[region]);
    }

    part.dpt = std::min(trkRes, caloRes);
    part.dphi = trkAngularResolution(part.pt, CHMS_[region], CHangular_[region]);
    part.deta = part.dphi;
}

void StandardParticleUncertainty::addUncertainty(particle& part, const jet& j){
    if(part.pdgid==211 || part.pdgid==11 || part.pdgid==13){ //charged particles
        addUncertaintyToCharged(part);
    } else if (part.pdgid==22){//photons
        addUncertaintyToNeutralEM(part);
    } else {
        addUncertaintyToNeutralHadron(part);
        if(part.pdgid!=130){
            std::cout << "Warning: unexpected pdgid " << part.pdgid << std::endl;
        }
    }
}

void StandardParticleUncertaintySmearedTracks::addUncertainty(particle& part, const jet& j){
    if(part.pdgid==13){ //charged particles
        addUncertaintyToCharged(part);
    } else if (part.pdgid==22 || part.pdgid==13){//photons + electrons
        addUncertaintyToNeutralEM(part);
    } else {
        addUncertaintyToNeutralHadron(part);
        if(part.pdgid!=130 && part.pdgid!=211){
            std::cout << "Warning: unexpected pdgid " << part.pdgid << std::endl;
        }
    }
}
