#ifndef SROTHMAN_MATCHING_V2_DELTARLIMITER_H
#define SROTHMAN_MATCHING_V2_DELTARLIMITER_H

#include "SRothman/SimonTools/src/jet.h"
#include <string>
#include <memory>

namespace matching {
    class DeltaRLimiter;
    using DeltaRLimiterPtr = std::unique_ptr<const DeltaRLimiter>;

    class DeltaRLimiter {
    public:
        virtual double evaluate(
                const double pt,
                const double eta,
                const double phi) const = 0;

        static DeltaRLimiterPtr get_deltaRlimiter(
                const std::string& mode,
                const double param1,
                const double param2);
    };

};

#endif
