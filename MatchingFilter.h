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

        static std::shared_ptr<MatchingFilter> getFilter(
                const std::string& behavior);
        static std::shared_ptr<MatchingFilter> getDRFilter(
                const std::string& behavior, 
                const std::vector<double>& etaBoundaries,
                const std::vector<double>& thresholds);
        static std::shared_ptr<MatchingFilter> getDRFilter(
                const std::string& behavior, 
                const std::vector<double>& etaBoundaries,
                const std::vector<double>& hardTerm, 
                const std::vector<double>& invTerm,
                const std::vector<double>& capTerm);
        static std::shared_ptr<MatchingFilter> getThresholdFilter(
                const std::vector<double>& boundaries, 
                const std::vector<double>& thresholds);
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

                const std::vector<std::string>& dRbehaviors,

                const std::vector<double>& EM0constDR,
                const std::vector<double>& EM0floatDR,
                const std::vector<double>& EM0capDR,

                const std::vector<double>& HAD0constDR,
                const std::vector<double>& HAD0floatDR,
                const std::vector<double>& HAD0capDR,

                const std::vector<double>& HADCHconstDR,
                const std::vector<double>& HADCHfloatDR,
                const std::vector<double>& HADCHcapDR,

                const std::vector<double>& ELEconstDR,
                const std::vector<double>& ELEfloatDR,
                const std::vector<double>& ELEcapDR,

                const std::vector<double>& MUconstDR,
                const std::vector<double>& MUfloatDR,
                const std::vector<double>& MUcapDR,

                const std::vector<double>& ECALEtaBoundaries,
                const std::vector<double>& HCALEtaBoundaries,
                const std::vector<double>& trkEtaBoundaries);

        bool pass(const particle& reco, const particle& gen);

        bool operator()(const particle& reco, const particle& gen){
            return pass(reco, gen);
        }

        static int pdgidToIndex(const particle& p);

    private:
        std::vector<std::shared_ptr<MatchingFilter>> softflavorfilters_;
        std::vector<std::shared_ptr<MatchingFilter>> hardflavorfilters_;
        std::vector<double> behaviorthresholds_;

        std::vector<std::shared_ptr<MatchingFilter>> chargefilters_;

        std::vector<std::shared_ptr<MatchingFilter>> dRfilters_;
        std::vector<std::shared_ptr<MatchingFilter>> threshfilters_;

};

#endif
