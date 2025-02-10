#include "refinePrefit.h"
#include "matchingUtil.h"
#include "SRothman/SimonTools/src/util.h"
#include "SRothman/SimonTools/src/isID.h"

class NONEprefitRefiner: public prefitRefiner{
    public:
        NONEprefitRefiner() {}
        ~NONEprefitRefiner() {}

        matchMap refine(const matchMap& recoToGen,
                        [[maybe_unused]] const simon::jet& recojet,
                        [[maybe_unused]] const simon::jet& genjet) const override {
            return invertMap(recoToGen);
        }
};

class OneGenOneRecoprefitRefiner: public prefitRefiner{
    public:
        OneGenOneRecoprefitRefiner() {}
        ~OneGenOneRecoprefitRefiner() {}

        matchMap refine(const matchMap& recoToGen,
                        const simon::jet& recojet,
                        const simon::jet& genjet) const override {
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

class OneGenOneRecoPerTypeprefitRefiner: public prefitRefiner{
    public:
        OneGenOneRecoPerTypeprefitRefiner() {}
        ~OneGenOneRecoPerTypeprefitRefiner() {}

        matchMap refine(const matchMap& recoToGen,
                        const simon::jet& recojet,
                        const simon::jet& genjet) const override {
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
                    if(isEM0(recopart)){
                        if(thischisq < bestChisqEM0){
                            bestChisqEM0 = thischisq;
                            bestIdxEM0 = iReco;
                        }
                    } else if(isHAD0(recopart)){
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

std::shared_ptr<prefitRefiner> prefitRefiner::get(const std::string& behavior){
    if(behavior == "None"){
        return std::make_shared<NONEprefitRefiner>();
    } else if(behavior == "OneGenOneReco"){
        return std::make_shared<OneGenOneRecoprefitRefiner>();
    } else if(behavior == "OneGenOneRecoPerType"){
        return std::make_shared<OneGenOneRecoPerTypeprefitRefiner>();
    } else {
        throw std::runtime_error("Unknown prefitRefinerType");
    }
}
