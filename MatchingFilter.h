#ifndef MATCHING_MATCHINGFILTER_H
#define MATCHING_MATCHINGFILTER_H

#include "SRothman/SimonTools/src/jets.h"
#include <memory>
#include <string> 

class MatchingFilter{
    public:
        MatchingFilter() {};
        virtual ~MatchingFilter(){};
        virtual bool pass(const particle& reco, const particle& gen) = 0;

        static std::shared_ptr<MatchingFilter> getFilter(const std::string& behavior);
};

class MatchingFilterEnsemble : public MatchingFilter{
    public:
        MatchingFilterEnsemble(
                const std::vector<std::string>& softflavorbehaviors,
                const std::vector<std::string>& hardflavorbehaviors,
                const std::vector<double>& behaviorthresholds,

                const std::vector<std::string>& chargebehaviors,

                const std::vector<double>& EM0thresholds,
                const std::vector<double>& HAD0thresholds,
                const std::vector<double>& HADCHthresholds,
                const std::vector<double>& ELEthresholds,
                const std::vector<double>& MUthresholds,

                const std::vector<double>& EM0dRcuts,
                const std::vector<double>& HAD0dRcuts,
                const std::vector<double>& HADCHdRcuts,
                const std::vector<double>& ELEdRcuts,
                const std::vector<double>& MUdRcuts,

                const std::vector<double>& ECALEtaBoundaries,
                const std::vector<double>& HCALEtaBoundaries,
                const std::vector<double>& trkEtaBoundaries);

        bool pass(const particle& reco, const particle& gen);

        bool operator()(const particle& reco, const particle& gen){
            return pass(reco, gen);
        }

        std::vector<std::shared_ptr<MatchingFilter>> softflavorfilters;
        std::vector<std::shared_ptr<MatchingFilter>> hardflavorfilters;
        std::vector<double> behaviorthresholds;

        std::vector<std::shared_ptr<MatchingFilter>> chargefilters;

        std::vector<std::shared_ptr<MatchingFilter>> dRfilters;
        std::vector<std::shared_ptr<MatchingFilter>> threshfilters;

        static int pdgidToIndex(int pdgid, int charge);
};

#endif
