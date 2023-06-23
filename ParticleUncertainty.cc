#include "ParticleUncertainty.h"
#include "SRothman/SimonTools/src/util.h"


void NaiveParticleUncertainty::addUncertainty(particle& part, const jet& j){
    part.dpt = 0.1 * part.pt;
    part.dphi = 0.05;
    part.deta = 0.05;
}

void RealisticParticleUncertainty::addUncertainty(particle& part, const jet& j){
    if(part.pdgid==211 || part.pdgid==11 || part.pdgid==13){ //charged particles
        if((part.pt < hardPt_) || (part.pdgid == 13)){ //muons get exemption
            addUncertaintyToCharged(part);
        } else {
            if(part.pdgid == 11){
                addUncertaintyToNeutralEM(part);
            } else {
                addUncertaintyToNeutralHadron(part);
            }
        }
    } else if (part.pdgid==22){//photons
        addUncertaintyToNeutralEM(part);
    } else if(part.pdgid==130){//neutral hadrons
        addUncertaintyToNeutralHadron(part);
    } else {
        std::cout << "Warning: unexpected pdgid " << part.pdgid << std::endl;
        addUncertaintyToNeutralHadron(part);
    }
}

void RealisticParticleUncertainty::addUncertaintyToNeutralEM(particle& part){
    int region = getEtaRegion(part.eta, ECALEtaBoundaries_);

    double E = part.pt * std::cosh(part.eta);
    part.dpt = std::sqrt( square(EMstochastic_[region]/std::sqrt(E)) 
                        + square(EMnoise_[region]/E)
                        + square(EMconstant_[region])) * part.pt;
    part.dphi = ECALgranularity_[region];
    part.deta = part.dphi;
}

void RealisticParticleUncertainty::addUncertaintyToNeutralHadron(particle& part){
    int region = getEtaRegion(part.eta, HCALEtaBoundaries_);

    double E = part.pt * std::cosh(part.eta);
    part.dpt = std::sqrt( square(HADstochastic_[region]/std::sqrt(E)) 
                        + square(HADconstant_[region])) * part.pt;
    part.dphi = HCALgranularity_[region];
    part.deta = part.dphi;
}

void RealisticParticleUncertainty::addUncertaintyToCharged(particle & part){
    int region = getEtaRegion(part.eta, trkEtaBoundaries_);

    part.dpt = std::sqrt( square(CHlinear_[region] * part.pt) 
                        + square(CHconstant_[region])) * part.pt;
    part.dphi = std::sqrt( square(CHMS_[region]/part.pt) 
                         + square(CHangular_[region]));
    part.deta = part.dphi;
}

int RealisticParticleUncertainty::getEtaRegion(double eta, const std::vector<double>& boundaries){
    for(unsigned i=0; i<boundaries.size(); ++i){
        if(std::abs(eta) < boundaries[i]){
            return i-1;
        }
    }
    return boundaries.size();
}
