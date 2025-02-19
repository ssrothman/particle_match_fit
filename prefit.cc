#include "prefit.h"
#include "SRothman/SimonTools/src/util.h"
#include "matchingUtil.h"
#include "SRothman/SimonTools/src/isID.h"

//#define DEBUG

class Prefitter{
    public:
        Prefitter(std::shared_ptr<MatchingFilter> filter) :
            filter_(filter) {}
        virtual ~Prefitter() {};

        std::vector<unsigned> operator()(const simon::particle& part, const simon::jet& j){
            return prefit(part, j);
        }

        virtual std::vector<unsigned> prefit(const simon::particle& part, const simon::jet& j)=0;

        static std::shared_ptr<Prefitter> get(const std::string& behavior, 
                            std::shared_ptr<MatchingFilter> filter);


    protected:
        std::shared_ptr<MatchingFilter> filter_;
};


class NOMATCHprefitter: public Prefitter{
    public:
        NOMATCHprefitter(std::shared_ptr<MatchingFilter> filter) :
            Prefitter(filter) {};
        ~NOMATCHprefitter() {};

        std::vector<unsigned> prefit([[maybe_unused]] const simon::particle& part, [[maybe_unused]] const simon::jet& j) override {
            return {};
        }
};

class CLOSESTprefitter : public Prefitter{
    public:
        CLOSESTprefitter(std::shared_ptr<MatchingFilter> filter):
            Prefitter(filter) {};
        ~CLOSESTprefitter() {};

        std::vector<unsigned> prefit(const simon::particle& part, const simon::jet& genjet) override {
            double minDR =std::numeric_limits<double>::infinity();
            int bestIdx = -1;

            for(unsigned i=0; i<genjet.nPart; ++i){
                const simon::particle& genpart = genjet.particles[i];
                if(filter_->pass(part, genpart)){
                    double dR = simon::deltaR2(part.eta, part.phi, 
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

class HARDESTprefitter : public Prefitter{
    public:
        HARDESTprefitter(std::shared_ptr<MatchingFilter> filter):
            Prefitter(filter) {};
        ~HARDESTprefitter() {};

        std::vector<unsigned> prefit(const simon::particle& part, const simon::jet& genjet) override {
            double maxPt = -std::numeric_limits<double>::infinity();
            int bestIdx = -1;

            for(unsigned i=0; i<genjet.nPart; ++i){
                const simon::particle& genpart = genjet.particles[i];
                if(filter_->pass(part, genpart)){
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

class BESTprefitter : public Prefitter{
    public:
        BESTprefitter(std::shared_ptr<MatchingFilter> filter):
            Prefitter(filter) {};
        ~BESTprefitter() {};

        std::vector<unsigned> prefit(const simon::particle& part, const simon::jet& j) override {
            double minChisq = std::numeric_limits<double>::infinity();
            int bestIdx = -1;

            for(unsigned i=0; i<j.nPart; ++i){
                const simon::particle& genpart = j.particles[i];
            

                if(filter_->pass(part, genpart)){
                    double chisq = chisquared(part, genpart, true);

                    if(chisq < minChisq || bestIdx==-1){
                        minChisq = chisq;
                        bestIdx = i;
                    }
#ifdef DEBUG
                    printf("\tgen %u: pass\n", i);
#endif
                } 
                else {
#ifdef DEBUG
                    printf("\tgen %u: fail\n", i);
#endif
                }
            }
            if(bestIdx<0){
                return {};
            } else {
                return {static_cast<unsigned>(bestIdx)};
            }
        }
};

class FLOATprefitter : public Prefitter{
    public:
        FLOATprefitter(std::shared_ptr<MatchingFilter> filter):
            Prefitter(filter) {};
        ~FLOATprefitter() {};

        std::vector<unsigned> prefit(const simon::particle& part, const simon::jet& j) override {
            std::vector<unsigned> result;
            for(unsigned i=0; i<j.nPart; ++i){
                const simon::particle& genpart = j.particles[i];
                if(filter_->pass(part, genpart)){
#ifdef DEBUG
                    printf("\tgen %u: pass\n", i);
#endif
                    result.push_back(i);
                } else {
#ifdef DEBUG
                    printf("\tgen %u: fail\n", i);
#endif
                }
            }
            return result;
        }
};

std::shared_ptr<Prefitter> Prefitter::get(const std::string& behavior, 
                                                    std::shared_ptr<MatchingFilter> filter){
    if(behavior == "NoMatch"){
        return std::make_shared<NOMATCHprefitter>(filter);
    } else if(behavior == "Closest"){
        return std::make_shared<CLOSESTprefitter>(filter);
    } else if(behavior == "Hardest"){
        return std::make_shared<HARDESTprefitter>(filter);
    } else if(behavior == "Best"){
        return std::make_shared<BESTprefitter>(filter);
    } else if(behavior == "Float"){
        return std::make_shared<FLOATprefitter>(filter);
    } else {
        throw std::invalid_argument("Invalid prefitter type");
    }
}

PrefitterEnsemble::PrefitterEnsemble(const std::vector<std::string>& behaviors, 
                                    std::shared_ptr<MatchingFilterEnsemble> filter){

    if(behaviors.size() != 5){
        throw std::invalid_argument("Must have inputs shape [EM0, HAD0, HADCH, ELE, MU]");
    }

    for(unsigned i=0; i<behaviors.size(); ++i){
        prefitters_.push_back(Prefitter::get(behaviors[i], filter));
    }
}

std::vector<unsigned> PrefitterEnsemble::prefit(const simon::particle& part, const simon::jet& j){
    if(isEM0(part)){
        return prefitters_[0]->prefit(part, j);
    } else if(isHAD0(part)){
        return prefitters_[1]->prefit(part, j);
    } else if(isHADCH(part)) {
        return prefitters_[2]->prefit(part, j);
    } else if(isELE(part)) {
        return prefitters_[3]->prefit(part, j);
    } else if(isMU(part)) {
        return prefitters_[4]->prefit(part, j);
    } else {
        printf("the bad pdgid, charge = %d, %d\n", part.pdgid, part.charge);
        throw std::invalid_argument("Invalid particle type in PrefitterEnsemble");
    }
}
