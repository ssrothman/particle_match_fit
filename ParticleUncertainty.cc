#include "ParticleUncertainty.h"
#include "SRothman/SimonTools/src/isID.h"
#include "SRothman/SimonTools/src/util.h"
#include "SRothman/SimonTools/src/etaRegion.h"

class NaiveParticleUncertainty : public ParticleUncertainty {
    public:
        NaiveParticleUncertainty() : ParticleUncertainty() {}
        ~NaiveParticleUncertainty() override {}
        void addUncertainty(simon::particle& part, const simon::jet& j) override;
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
        void addUncertaintyToNeutralEM(simon::particle& part);
        void addUncertaintyToNeutralHadron(simon::particle& part);
        void addUncertaintyToCharged(simon::particle& part);

        std::vector<double> EMstochastic_, EMconstant_;
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
                    EMstochastic, EMconstant,
                    ECALgranularity, ECALEtaBoundaries, 
                    HADstochastic, HADconstant, 
                    HCALgranularity, HCALEtaBoundaries,
                    CHlinear, CHconstant, 
                    CHMS, CHangular,
                    trkEtaBoundaries) {}
        ~StandardParticleUncertainty() override {}
        void addUncertainty(simon::particle& part, const simon::jet& j) override;
};

class StandardParticleUncertaintySmearedTracks : public RealisticParticleUncertainty {
    public:
        StandardParticleUncertaintySmearedTracks(
                const std::vector<double>& EMstochastic,
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
                    EMstochastic, EMconstant,
                    ECALgranularity, ECALEtaBoundaries, 
                    HADstochastic, HADconstant, 
                    HCALgranularity, HCALEtaBoundaries,
                    CHlinear, CHconstant, 
                    CHMS, CHangular,
                    trkEtaBoundaries) {}
        ~StandardParticleUncertaintySmearedTracks() override {}
        void addUncertainty(simon::particle& part, const simon::jet& j) override;
};

static double caloResolution(const double& eta, const double& pt, const double& A, const double& C){
    double E = pt * std::cosh(eta);
    return std::sqrt(simon::square(A/std::sqrt(E)) + simon::square(C)) * pt;
}

static double trkResolution(const double& pt, const double& A, const double& B){
    return std::sqrt(simon::square(A*pt) + simon::square(B)) * pt;
}

static double trkAngularResolution(const double& pt, const double& A, const double& B){
    return std::sqrt(simon::square(A/pt) + simon::square(B));
}

void NaiveParticleUncertainty::addUncertainty(simon::particle& part, [[maybe_unused]] const simon::jet& j){
    part.dpt = 0.1 * part.pt;
    part.dphi = 0.05;
    part.deta = 0.05;
}

void RealisticParticleUncertainty::addUncertaintyToNeutralEM(simon::particle& part){
    int region = simon::getEtaRegion(part.eta, ECALEtaBoundaries_);

    part.dpt = caloResolution(part.eta, part.pt, EMstochastic_[region], EMconstant_[region]);
    part.dphi = ECALgranularity_[region];
    part.deta = part.dphi;
}

void RealisticParticleUncertainty::addUncertaintyToNeutralHadron(simon::particle& part){
    int region = simon::getEtaRegion(part.eta, HCALEtaBoundaries_);

    part.dpt = caloResolution(part.eta, part.pt, HADstochastic_[region], HADconstant_[region]);
    part.dphi = HCALgranularity_[region];
    part.deta = part.dphi;
}

void RealisticParticleUncertainty::addUncertaintyToCharged(simon::particle& part){
    int region = simon::getEtaRegion(part.eta, trkEtaBoundaries_);

    double trkRes = trkResolution(part.pt, CHlinear_[region], CHconstant_[region]);
    double caloRes = 9e9;
    if(isELE(part)){
        caloRes = caloResolution(part.eta, part.pt, 
                                        EMstochastic_[region], 
                                        EMconstant_[region]);
    } else if(isMU(part)){//muons are tracker-only objects
        caloRes = caloResolution(part.eta, part.pt, 
                                HADstochastic_[region], 
                                HADconstant_[region]);
    }

    part.dpt = std::min(trkRes, caloRes);
    part.dphi = trkAngularResolution(part.pt, CHMS_[region], CHangular_[region]);
    part.deta = part.dphi;
}

void StandardParticleUncertainty::addUncertainty(simon::particle& part, [[maybe_unused]] const simon::jet& j){
    if(part.charge!=0){ //charged particles
        addUncertaintyToCharged(part);
    } else if (isEM0(part)){//photons
        addUncertaintyToNeutralEM(part);
    } else {
        addUncertaintyToNeutralHadron(part);
    }
}

void StandardParticleUncertaintySmearedTracks::addUncertainty(simon::particle& part, [[maybe_unused]] const simon::jet& j){
    if(isMU(part)){ //muons can't be smeared
        addUncertaintyToCharged(part);
    } else if (isEM0(part) || isELE(part)){//photons + electrons
        addUncertaintyToNeutralEM(part);
    } else {
        addUncertaintyToNeutralHadron(part);
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
        const std::vector<double>& EMconstant,
        const std::vector<double>& ECALgranularityEta,
        [[maybe_unused]] const std::vector<double>& ECALgranularityPhi,
        const std::vector<double>& ECALEtaBoundaries,
        const std::vector<double>& HADstochastic,
        const std::vector<double>& HADconstant,
        const std::vector<double>& HCALgranularityEta,
        [[maybe_unused]] const std::vector<double>& HCALgranularityPhi,
        const std::vector<double>& HCALEtaBoundaries,
        const std::vector<double>& CHlinear,
        const std::vector<double>& CHconstant,
        const std::vector<double>& CHMSeta,
        [[maybe_unused]] const std::vector<double>& CHMSphi,
        const std::vector<double>& CHangularEta,
        [[maybe_unused]] const std::vector<double>& CHangularPhi,
        const std::vector<double>& trkEtaBoundaries){

    if(behavior == "Naive"){
        return std::make_shared<NaiveParticleUncertainty>();
    } else if(behavior == "Standard"){
        return std::make_shared<StandardParticleUncertainty>(
                EMstochastic, EMconstant,
                ECALgranularityEta, ECALEtaBoundaries, 
                HADstochastic, HADconstant, 
                HCALgranularityEta, HCALEtaBoundaries,
                CHlinear, CHconstant, 
                CHMSeta, CHangularEta,
                trkEtaBoundaries);
    } else if(behavior == "SmearedTracks"){
        return std::make_shared<StandardParticleUncertaintySmearedTracks>(
                EMstochastic, EMconstant,
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
