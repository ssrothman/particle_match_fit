#include "prefit.h"
#include "SRothman/SimonTools/src/util.h"
#include "matchingUtil.h"

class NOMATCHprefitter: public prefitter{
    public:
        NOMATCHprefitter(std::shared_ptr<MatchingFilter> filter) :
            prefitter(filter) {};
        ~NOMATCHprefitter() {};

        std::vector<unsigned> operator()(const particle& part, const jet& j) override {
            return {};
        }
};

class CLOSESTprefitter : public prefitter{
    public:
        CLOSESTprefitter(std::shared_ptr<MatchingFilter> filter):
            prefitter(filter) {};
        ~CLOSESTprefitter() {};

        std::vector<unsigned> operator()(const particle& part, const jet& genjet) override {
            double minDR = 99999;
            int bestIdx = -1;

            for(unsigned i=0; i<genjet.nPart; ++i){
                const particle& genpart = genjet.particles[i];
                if(filter_->allowMatch(part, genpart, genjet)){
                    double dR = dR2(part.eta, part.phi, 
                                    genpart.eta, genpart.phi);
                    if(dR<minDR || bestIdx==-1){
                        minDR = dR;
                        bestIdx = i;
                    }
                } 
            }
            if(bestIdx<0){
                return {};
            } else {
                return {static_cast<unsigned>(bestIdx)};
            }
        }
};

class HARDESTprefitter : public prefitter{
    public:
        HARDESTprefitter(std::shared_ptr<MatchingFilter> filter):
            prefitter(filter) {};
        ~HARDESTprefitter() {};

        std::vector<unsigned> operator()(const particle& part, const jet& genjet) override {
            double maxPt = 0;
            int bestIdx = -1;

            for(unsigned i=0; i<genjet.nPart; ++i){
                const particle& genpart = genjet.particles[i];
                if(filter_->allowMatch(part, genpart, genjet)){
                    if(genpart.pt>maxPt || bestIdx==-1){
                        maxPt = genpart.pt;
                        bestIdx = i;
                    }
                } 
            }
            if(bestIdx<0){
                return {};
            } else {
                return {static_cast<unsigned>(bestIdx)};
            }
        }
};

class BESTprefitter : public prefitter{
    public:
        BESTprefitter(std::shared_ptr<MatchingFilter> filter):
            prefitter(filter) {};
        ~BESTprefitter() {};

        std::vector<unsigned> operator()(const particle& part, const jet& j) override {
            double minChisq = std::numeric_limits<double>::infinity();
            int bestIdx = -1;

            for(unsigned i=0; i<j.nPart; ++i){
                const particle& genpart = j.particles[i];
            

                if(filter_->allowMatch(part, genpart, j)){
                    double chisq = chisquared(part, genpart, true);

                    if(chisq < minChisq || bestIdx==-1){
                        minChisq = chisq;
                        bestIdx = i;
                    }
                } 
            }
            if(bestIdx<0){
                return {};
            } else {
                return {static_cast<unsigned>(bestIdx)};
            }
        }
};

class FLOATprefitter : public prefitter{
    public:
        FLOATprefitter(std::shared_ptr<MatchingFilter> filter):
            prefitter(filter) {};
        ~FLOATprefitter() {};

        std::vector<unsigned> operator()(const particle& part, const jet& j) override {
            std::vector<unsigned> result;
            for(unsigned i=0; i<j.nPart; ++i){
                const particle& genpart = j.particles[i];
                if(filter_->allowMatch(part, genpart, j)){
                    result.push_back(i);
                } 
            }
            return result;
        }
};

std::shared_ptr<prefitter> prefitter::getPrefitter(const enum prefitterType& behavior, 
                                                    std::shared_ptr<MatchingFilter> filter){
    switch(behavior){
        case prefitterType::NOMATCH:
            return std::make_shared<NOMATCHprefitter>(filter);
        case prefitterType::CLOSEST:
            return std::make_shared<CLOSESTprefitter>(filter);
        case prefitterType::HARDEST:
            return std::make_shared<HARDESTprefitter>(filter);
        case prefitterType::BEST:
            return std::make_shared<BESTprefitter>(filter);
        case prefitterType::FLOAT:
            return std::make_shared<FLOATprefitter>(filter);
        default:
            throw std::invalid_argument("Invalid prefitter type");
    }
}

