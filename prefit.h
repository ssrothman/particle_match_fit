#ifndef PREFIT_H
#define PREFIT_H

#include <vector>
#include "SRothman/SimonTools/src/deltaR.h"
#include "MatchingFilter.h"
#include "SRothman/SimonTools/src/jets.h"
#include <memory>
#include <string>

class Prefitter;

class PrefitterEnsemble{
    public:
        PrefitterEnsemble(const std::vector<std::string>& behaviors, std::shared_ptr<MatchingFilterEnsemble> filters);
        
        std::vector<unsigned> prefit(const particle& part, const jet& j);

        std::vector<unsigned> operator()(const particle& part, const jet& j){
            return prefit(part, j);
        }


    private:
        std::vector<std::shared_ptr<Prefitter>> prefitters_;
};

#endif
