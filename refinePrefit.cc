#include "refinePrefit.h"
#include "matchingUtil.h"
#include "SRothman/SimonTools/src/util.h"

class NONEprefitRefiner: public prefitRefiner{
    public:
        NONEprefitRefiner() {}
        ~NONEprefitRefiner() {}

        matchMap refine(const matchMap& recoToGen,
                        const jet& recojet,
                        const jet& genjet) const override {
            return invertMap(recoToGen);
        }
};

class GLOBAL_SINGLE_GENprefitRefiner: public prefitRefiner{
    public:
        GLOBAL_SINGLE_GENprefitRefiner() {}
        ~GLOBAL_SINGLE_GENprefitRefiner() {}

        matchMap refine(const matchMap& recoToGen,
                        const jet& recojet,
                        const jet& genjet) const override {
            matchMap genToReco = invertMap(recoToGen);

            matchMap result;

            for(auto& matches : genToReco){
                double bestChisq = std::numeric_limits<double>::infinity();
                int bestIdx = -1;

                const auto& genpart = genjet.particles[matches.first];
                for(unsigned iReco : matches.second){
                    const auto& recopart = recojet.particles[iReco];
                    double thischisq = chisquared(recopart,genpart, false);
                    if(thischisq<bestChisq){
                        bestChisq = thischisq;
                        bestIdx = iReco;
                    }
                }

                if(bestIdx>=0)
                    result[matches.first] = {static_cast<unsigned>(bestIdx)};
                else
                    result[matches.first] = {};
            }

            return result;
        }
};

class SINGLE_GENprefitRefiner: public prefitRefiner{
    public:
        SINGLE_GENprefitRefiner() {}
        ~SINGLE_GENprefitRefiner() {}

        matchMap refine(const matchMap& recoToGen,
                        const jet& recojet,
                        const jet& genjet) const override {
            matchMap genToReco = invertMap(recoToGen);

            matchMap result;

            for(auto& matches : genToReco){
                double bestChisqCH = std::numeric_limits<double>::infinity();
                double bestChisqEM0 = std::numeric_limits<double>::infinity();
                double bestChisqHAD0 = std::numeric_limits<double>::infinity();

                int bestIdxCH = -1;
                int bestIdxEM0 = -1;
                int bestIdxHAD0 = -1;

                const auto& genpart = genjet.particles[matches.first];
                for(unsigned iReco : matches.second){
                    const auto& recopart = recojet.particles[iReco];
                    double thischisq = chisquared(recopart,genpart, false);
                    if(recopart.pdgid == 22){
                        if(thischisq < bestChisqEM0){
                            bestChisqEM0 = thischisq;
                            bestIdxEM0 = iReco;
                        }
                    } else if(recopart.pdgid == 130){
                        if(thischisq < bestChisqHAD0){
                            bestChisqHAD0 = thischisq;
                            bestIdxHAD0 = iReco;
                        }
                    }
                    else if(recopart.charge != 0){
                        if(thischisq < bestChisqCH){
                            bestChisqCH = thischisq;
                            bestIdxCH = iReco;
                        }
                    } else {
                        throw std::runtime_error("unknown particle type");
                    }
                }

                std::vector<unsigned> bestIdxs;
                if(bestIdxCH>=0)
                    bestIdxs.push_back(bestIdxCH);
                if(bestIdxEM0>=0)
                    bestIdxs.push_back(bestIdxEM0);
                if(bestIdxHAD0>=0)
                    bestIdxs.push_back(bestIdxHAD0);
                
                result[matches.first] = bestIdxs;
            }

            return result;
        }
};

prefitRefiner::matchMap prefitRefiner::invertMap(const matchMap& m) const {
    matchMap inverted;
    for(const auto& pair : m){
        for(unsigned i : pair.second){
            inverted[i].push_back(pair.first);
        }
    }
    return inverted;
}

std::shared_ptr<prefitRefiner> prefitRefiner::getRefiner(const enum prefitRefinerType& behavior){
    switch(behavior){
        case prefitRefinerType::NONE:
            return std::make_shared<NONEprefitRefiner>();
        case prefitRefinerType::GLOBAL_SINGLE_GEN:
            return std::make_shared<GLOBAL_SINGLE_GENprefitRefiner>();
        case prefitRefinerType::SINGLE_GEN:
            return std::make_shared<SINGLE_GENprefitRefiner>();
        default:
            throw std::runtime_error("Unknown prefitRefinerType");
    }
}
