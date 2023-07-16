#ifndef MATCHING_PARTICLEUNCERTAINTY_H
#define MATCHING_PARTICLEUNCERTAINTY_H

#include "SRothman/SimonTools/src/jets.h"
#include <vector>
#include <memory>

enum class uncertaintyType{
    NAIVE = 0,
    STANDARD = 1,
    SMEAREDTRACKS = 2
};

class ParticleUncertainty{
    public:
        ParticleUncertainty(){};
        virtual ~ParticleUncertainty(){};
        virtual void addUncertainty(particle& part, const jet& j) = 0;

    static std::shared_ptr<ParticleUncertainty> getUncertainty(
            const enum uncertaintyType& behavior);

    static std::shared_ptr<ParticleUncertainty> getUncertainty(
            const enum uncertaintyType& behavior, 
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
            const std::vector<double>& trkEtaBoundaries);
};

#endif
