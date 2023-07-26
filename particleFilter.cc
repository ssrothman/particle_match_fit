#include "particleFilter.h"

class ALLparticleFilter : public particleFilter{
public:
    ALLparticleFilter() {};
    ~ALLparticleFilter() {};
    bool pass(const particle& part) override{
        return true;
    }
};

class NONEparticleFilter : public particleFilter{
public:
    NONEparticleFilter() {};
    ~NONEparticleFilter() {};
    bool pass(const particle& part) override{
        return false;
    }
};

class CHARGEDparticleFilter : public particleFilter{
public:
    CHARGEDparticleFilter() {};
    ~CHARGEDparticleFilter() {};
    bool pass(const particle& part) override{
        return part.charge != 0;
    }
};

class NEUTRALparticleFilter : public particleFilter{
public:
    NEUTRALparticleFilter() {};
    ~NEUTRALparticleFilter() {};
    bool pass(const particle& part) override{
        return part.charge == 0;
    }
};

class EM0particleFilter : public particleFilter{
public:
    EM0particleFilter() {};
    ~EM0particleFilter() {};
    bool pass(const particle& part) override{
        return part.pdgid == 22;
    }
};

class HAD0particleFilter : public particleFilter{
public:
    HAD0particleFilter() {};
    ~HAD0particleFilter() {};
    bool pass(const particle& part) override{
        return part.pdgid == 130;
    }
};

std::shared_ptr<particleFilter> particleFilter::getParticleFilter(particleFilterType type){
    switch(type){
        case particleFilterType::ALL:
            return std::make_shared<ALLparticleFilter>();
        case particleFilterType::NONE:
            return std::make_shared<NONEparticleFilter>();
        case particleFilterType::CHARGED:
            return std::make_shared<CHARGEDparticleFilter>();
        case particleFilterType::NEUTRAL:
            return std::make_shared<NEUTRALparticleFilter>();
        case particleFilterType::EM0:
            return std::make_shared<EM0particleFilter>();
        case particleFilterType::HAD0:
            return std::make_shared<HAD0particleFilter>();
        default:
            throw std::runtime_error("particleFilter::getParticleFilter: Unknown particleFilterType");
    }
}
