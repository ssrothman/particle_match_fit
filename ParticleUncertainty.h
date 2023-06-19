#ifndef MATCHING_PARTICLEUNCERTAINTY_H
#define MATCHING_PARTICLEUNCERTAINTY_H

#include "SRothman/Matching/src/toyjets/common.h"
#include <vector>

class ParticleUncertainty{
    public:
        ParticleUncertainty(){};
        virtual ~ParticleUncertainty(){};
        virtual void addUncertainty(particle& part, const jet& j) = 0;
};

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
     *
     * NB resolutions are treated to be constant as a function of eta within the tracker volume
     *  This is a pretty bad assumption, and might need to be updated
     *
     * Additionally, we take into account poor tracking performance at high pT, which is particularly seen in jets
     * if the particle has pt > hardPt_ then we revert to calorimeter resolutions 
     * NB muons are excepted from this, as they are not calorimeter objects at all
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
                const std::vector<double>& trkEtaBoundaries,
                double hardPt
                ) :
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
            trkEtaBoundaries_(trkEtaBoundaries),
            hardPt_(hardPt) {}
        ~RealisticParticleUncertainty() override {};
        void addUncertainty(particle& part, const jet& j) override;
    private:
        void addUncertaintyToNeutralEM(particle& part);
        void addUncertaintyToNeutralHadron(particle& part);
        void addUncertaintyToCharged(particle& part);

        int getEtaRegion(double eta, const std::vector<double>& boundaries);

        std::vector<double> EMstochastic_, EMnoise_, EMconstant_;
        std::vector<double> ECALgranularity_;
        std::vector<double> ECALEtaBoundaries_;

        std::vector<double> HADstochastic_, HADconstant_;
        std::vector<double> HCALgranularity_;
        std::vector<double> HCALEtaBoundaries_;

        std::vector<double> CHlinear_, CHconstant_;
        std::vector<double> CHMS_, CHangular_;
        std::vector<double> trkEtaBoundaries_;

        double hardPt_;
};

#endif
