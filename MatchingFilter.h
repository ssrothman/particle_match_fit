#ifndef MATCHING_MATCHINGFILTER_H
#define MATCHING_MATCHINGFILTER_H

#include "SRothman/SimonTools/src/jets.h"
#include <memory>
#include <string>

class MatchingFilter{
    public:
        MatchingFilter(double threshold) : threshold_(threshold) {};
        virtual ~MatchingFilter(){};
        virtual bool pass(const particle& reco, const particle& gen) = 0;

        static std::shared_ptr<MatchingFilter> getFilter(const std::string& behavior, double threshold);

        double threshold_;
};

class MatchingFilterEnsemble{
    public:
        MatchingFilterEnsemble(std::vector<std::string> behaviors, std::vector<double> thresholds);
        
        bool pass(const particle& reco, const particle& gen);

        bool operator()(const particle& reco, const particle& gen){
            return pass(reco, gen);
        }

        std::vector<std::shared_ptr<MatchingFilter>> filters;
};

#endif
