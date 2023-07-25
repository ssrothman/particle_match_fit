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
        prefitter(std::shared_ptr<MatchingFilter> filter) :
            filter_(filter) {}
        virtual ~prefitter() {};

        virtual std::vector<unsigned> operator()(const particle& part, const jet& j) = 0;

        static std::shared_ptr<prefitter> getPrefitter(const enum prefitterType& behavior, 
                            std::shared_ptr<MatchingFilter> filter);


    protected:
        std::shared_ptr<MatchingFilter> filter_;
};

#endif
