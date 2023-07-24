#ifndef REFINE_PREFIT_H
#define REFINE_PREFIT_H

#include <vector>
#include <memory>
#include <armadillo>
#include "SRothman/SimonTools/src/jets.h"

enum class prefitRefinerType{
    NONE=0,
    GLOBAL_SINGLE_GEN=1,
    SINGLE_GEN=2,
};

class prefitRefiner{
    public:
        prefitRefiner() {}
        virtual ~prefitRefiner() {}

        typedef std::unordered_map<unsigned, std::vector<unsigned>> matchMap;

        virtual matchMap refine(const matchMap& recoToGen,
                                const jet& recojet,
                                const jet& genjet) const = 0;
                
        static std::shared_ptr<prefitRefiner> getRefiner(const enum prefitRefinerType& behavior);

    protected:
        matchMap invertMap(const matchMap& m) const;

};
#endif
