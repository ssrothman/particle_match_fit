#ifndef MATCHING_PARTICLEUNCERTAINTY_H
#define MATCHING_PARTICLEUNCERTAINTY_H

#include "SRothman/SimonTools/src/jet.h"
#include <vector>
#include <memory>
#include <string>

class ParticleUncertainty{
    public:
        ParticleUncertainty(){};
        virtual ~ParticleUncertainty(){};
        virtual void addUncertainty(simon::particle& part, const simon::jet& j) = 0;

    static std::shared_ptr<ParticleUncertainty> get(
            const std::string& behavior);

    static std::shared_ptr<ParticleUncertainty> get(
            const std::string& behavior, 
            const std::vector<double>& EMstochastic, 
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
            const std::vector<double>& trkEtaBoundaries);
};

#endif
