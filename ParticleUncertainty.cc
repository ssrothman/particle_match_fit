#include "ParticleUncertainty.h"
#include "SRothman/SimonTools/src/util.h"

class NaiveParticleUncertainty : public ParticleUncertainty {
    public:
        NaiveParticleUncertainty() : ParticleUncertainty() {}
        ~NaiveParticleUncertainty() override {}
        void addUncertainty(particle& part, const jet& j) override;
};


class RealisticParticleUncertainty : public ParticleUncertainty{
    /*
     * Attempt to emulate something realistic in CMS
     *
     * For neutral EM particles, resolutions are determined by ECAL
     *  sigma(E)/E = A/sqrt(E) + B/E + C
     *  sigma(phi) = sigma(eta) = D
     *
     * For neutral hadrons, resolutions are determined by HCAL
     *  sigma(E)/E = A/sqrt(E) + C
     *  sigma(phi) = sigma(eta) = D
     *
     * For charged particles, resolutions are determined by tracker
     *   sigma(pT)/pT = A*pT + B
     *   sigma(phi) = sigma(eta) = C/pT + D
     */
    public:
        RealisticParticleUncertainty(
                const std::vector<double>& EMstochastic,
                const std::vector<double>& EMnoise,
                const std::vector<double>& EMconstant,
                const std::vector<double>& ECALgranularity,
                const std::vector<double>& ECALEtaBoundaries,
                const std::vector<double>& HADstochastic,
                const std::vector<double>& HADconstant,
                const std::vector<double>& HCALgranularity,
                const std::vector<double>& HCALEtaBoundaries,
                const std::vector<double>& CHlinear,
                const std::vector<double>& CHconstant,
                const std::vector<double>& CHMS,
                const std::vector<double>& CHangular,
                const std::vector<double>& trkEtaBoundaries) :
            EMstochastic_(EMstochastic),
            EMnoise_(EMnoise),
            EMconstant_(EMconstant),
            ECALgranularity_(ECALgranularity),
            ECALEtaBoundaries_(ECALEtaBoundaries),
            HADstochastic_(HADstochastic),
            HADconstant_(HADconstant),
            HCALgranularity_(HCALgranularity),
            HCALEtaBoundaries_(HCALEtaBoundaries),
            CHlinear_(CHlinear),
            CHconstant_(CHconstant),
            CHMS_(CHMS),
            CHangular_(CHangular),
            trkEtaBoundaries_(trkEtaBoundaries) {}
    protected:
        void addUncertaintyToNeutralEM(particle& part);
        void addUncertaintyToNeutralHadron(particle& part);
        void addUncertaintyToCharged(particle& part);

        std::vector<double> EMstochastic_, EMnoise_, EMconstant_;
        std::vector<double> ECALgranularity_;
        std::vector<double> ECALEtaBoundaries_;

        std::vector<double> HADstochastic_, HADconstant_;
        std::vector<double> HCALgranularity_;
        std::vector<double> HCALEtaBoundaries_;

        std::vector<double> CHlinear_, CHconstant_;
        std::vector<double> CHMS_, CHangular_;
        std::vector<double> trkEtaBoundaries_;
};

class StandardParticleUncertainty : public RealisticParticleUncertainty {
    public:
        StandardParticleUncertainty(
                const std::vector<double>& EMstochastic,
                const std::vector<double>& EMnoise,
                const std::vector<double>& EMconstant,
                const std::vector<double>& ECALgranularity,
                const std::vector<double>& ECALEtaBoundaries,
                const std::vector<double>& HADstochastic,
                const std::vector<double>& HADconstant,
                const std::vector<double>& HCALgranularity,
                const std::vector<double>& HCALEtaBoundaries,
                const std::vector<double>& CHlinear,
                const std::vector<double>& CHconstant,
                const std::vector<double>& CHMS,
                const std::vector<double>& CHangular,
                const std::vector<double>& trkEtaBoundaries) :
            RealisticParticleUncertainty(
                    EMstochastic, EMnoise, EMconstant,
                    ECALgranularity, ECALEtaBoundaries, 
                    HADstochastic, HADconstant, 
                    HCALgranularity, HCALEtaBoundaries,
                    CHlinear, CHconstant, 
                    CHMS, CHangular,
                    trkEtaBoundaries) {}
        ~StandardParticleUncertainty() override {}
        void addUncertainty(particle& part, const jet& j) override;
};

class StandardParticleUncertaintySmearedTracks : public RealisticParticleUncertainty {
    public:
        StandardParticleUncertaintySmearedTracks(
                const std::vector<double>& EMstochastic,
                const std::vector<double>& EMnoise,
                const std::vector<double>& EMconstant,
                const std::vector<double>& ECALgranularity,
                const std::vector<double>& ECALEtaBoundaries,
                const std::vector<double>& HADstochastic,
                const std::vector<double>& HADconstant,
                const std::vector<double>& HCALgranularity,
                const std::vector<double>& HCALEtaBoundaries,
                const std::vector<double>& CHlinear,
                const std::vector<double>& CHconstant,
                const std::vector<double>& CHMS,
                const std::vector<double>& CHangular,
                const std::vector<double>& trkEtaBoundaries) :
            RealisticParticleUncertainty(
                    EMstochastic, EMnoise, EMconstant,
                    ECALgranularity, ECALEtaBoundaries, 
                    HADstochastic, HADconstant, 
                    HCALgranularity, HCALEtaBoundaries,
                    CHlinear, CHconstant, 
                    CHMS, CHangular,
                    trkEtaBoundaries) {}
        ~StandardParticleUncertaintySmearedTracks() override {}
        void addUncertainty(particle& part, const jet& j) override;
};



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

std::shared_ptr<ParticleUncertainty> ParticleUncertainty::get(const std::string& behavior){
    if(behavior == "Naive"){
        return std::make_shared<NaiveParticleUncertainty>();
    } else if(behavior == "Standard" || behavior == "SmearedTracks"){
        throw std::invalid_argument("StandardParticleUncertainty requires more parameters");
    } else {
        throw std::invalid_argument("Unknown uncertainty type");
    }
}


std::shared_ptr<ParticleUncertainty> ParticleUncertainty::get(
        const std::string& behavior,
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
        const std::vector<double>& trkEtaBoundaries){

    if(behavior == "Naive"){
        return std::make_shared<NaiveParticleUncertainty>();
    } else if(behavior == "Standard"){
        return std::make_shared<StandardParticleUncertainty>(
                EMstochastic, EMnoise, EMconstant,
                ECALgranularityEta, ECALEtaBoundaries, 
                HADstochastic, HADconstant, 
                HCALgranularityEta, HCALEtaBoundaries,
                CHlinear, CHconstant, 
                CHMSeta, CHangularEta,
                trkEtaBoundaries);
    } else if(behavior == "SmearedTracks"){
        return std::make_shared<StandardParticleUncertaintySmearedTracks>(
                EMstochastic, EMnoise, EMconstant,
                ECALgranularityEta, ECALEtaBoundaries, 
                HADstochastic, HADconstant, 
                HCALgranularityEta, HCALEtaBoundaries,
                CHlinear, CHconstant, 
                CHMSeta, CHangularEta,
                trkEtaBoundaries);
    } else {
        throw std::invalid_argument("Unknown uncertainty type");
    }
}
