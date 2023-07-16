#ifndef MATCHING_MATCHINGFILTER_H
#define MATCHING_MATCHINGFILTER_H

#include "SRothman/SimonTools/src/jets.h"
#include <memory>

enum class matchFilterType{
    DR = 0,
    CHARGESIGN = 1,
    CHARGE = 2,
    REALISTIC = 3,
    LOSTTRACK = 4
};

class MatchingFilter{
    public:
        MatchingFilter(double threshold) : threshold_(threshold) {};
        virtual ~MatchingFilter(){};
        virtual bool allowMatch(const particle& reco, const particle& gen, const jet& j) = 0;

        static std::shared_ptr<MatchingFilter> getFilter(const enum matchFilterType& behavior, double threshold);
        static std::shared_ptr<MatchingFilter> getFilter(const enum matchFilterType& behavior, double threshold, double softPt, double hardPt);

    protected:
        bool passDR(const particle& reco, const particle& gen, const jet& j);
        bool sameCharge(const particle& reco, const particle& gen, const jet& j);
        bool sameChargeSign(const particle& reco, const particle& gen, const jet& j);
    private:
        double threshold_;
};

#endif
