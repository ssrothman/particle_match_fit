#ifndef PREFIT_H
#define PREFIT_H

#include <vector>
#include "SRothman/SimonTools/src/deltaR.h"
#include "MatchingFilter.h"
#include "SRothman/SimonTools/src/jets.h"
#include <memory>

enum class prefitterType{
    NOMATCH = 0,
    CLOSEST = 1,
    HARDEST = 2,
    BEST = 3,
    FLOAT = 4,
};

class prefitter{
    public:
        prefitter(std::shared_ptr<MatchingFilter> filter, const std::vector<bool>& excludeGen) :
            filter_(filter),
            excludeGen_(excludeGen) {}
        virtual ~prefitter() {};

        virtual std::vector<unsigned> operator()(const particle& part, const jet& j) = 0;

        static std::shared_ptr<prefitter> getPrefitter(const enum prefitterType& behavior, 
                            std::shared_ptr<MatchingFilter> filter,
                            const std::vector<bool>& excludeGen);

    protected:
        std::shared_ptr<MatchingFilter> filter_;
        std::vector<bool> excludeGen_;
};

#endif
