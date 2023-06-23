#ifndef MATCHING_MATCHINGFILTER_H
#define MATCHING_MATCHINGFILTER_H

#include "SRothman/SimonTools/src/jets.h"

class MatchingFilter{
    public:
        MatchingFilter(double threshold) : threshold_(threshold) {};
        virtual ~MatchingFilter(){};
        virtual bool allowMatch(const particle& part1, const particle& part2, const jet& j) = 0;
    protected:
        bool passDR(const particle& part1, const particle& part2, const jet& j);
        bool sameCharge(const particle& part1, const particle& part2, const jet& j);
        bool sameChargeSign(const particle& part1, const particle& part2, const jet& j);
    private:
        double threshold_;
};

class DRFilter : public MatchingFilter {
    public:
        DRFilter(double threshold) : MatchingFilter(threshold) {}
        ~DRFilter() override {};
        bool allowMatch(const particle& part1, const particle& part2, const jet& j) override;
};

class ChargeSignFilter : public MatchingFilter {
    public:
        ChargeSignFilter(double threshold) : MatchingFilter(threshold) {};
        ~ChargeSignFilter() override {};
        bool allowMatch(const particle& part1, const particle& part2, const jet& j) override;
};

class ChargeFilter : public MatchingFilter {
    public:
        ChargeFilter(double threshold) : MatchingFilter(threshold) {};
        ~ChargeFilter() override {};
        bool allowMatch(const particle& part1, const particle& part2, const jet& j) override;
};

class RealisticFilter : public MatchingFilter {
    public:
        RealisticFilter(double threshold, double softPt, double hardPt) : 
            MatchingFilter(threshold),
            softPt_(softPt), hardPt_(hardPt) {};
        ~RealisticFilter() override {};
        bool allowMatch(const particle& part1, const particle& part2, const jet& j) override;
    private:
        double softPt_, hardPt_;
};

#endif
