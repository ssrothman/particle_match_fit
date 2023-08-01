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

class EM0HAD0particleFilter : public particleFilter{
public:
    EM0HAD0particleFilter() {};
    ~EM0HAD0particleFilter() {};
    bool pass(const particle& part) override{
        return part.pdgid == 22 || part.pdgid == 130;
    }
};

class EM0CHARGEDparticleFilter : public particleFilter{
public:
    EM0CHARGEDparticleFilter() {};
    ~EM0CHARGEDparticleFilter() {};
    bool pass(const particle& part) override{
        return part.pdgid == 22 || part.charge != 0;
    }
};

class HAD0CHARGEDparticleFilter : public particleFilter{
public:
    HAD0CHARGEDparticleFilter() {};
    ~HAD0CHARGEDparticleFilter() {};
    bool pass(const particle& part) override{
        return part.pdgid == 130 || part.charge != 0;
    }
};

std::shared_ptr<particleFilter> particleFilter::getParticleFilter(const std::string& type){
    if(type == "ALL"){
        return std::make_shared<ALLparticleFilter>();
    } else if(type == "NONE"){
        return std::make_shared<NONEparticleFilter>();
    } else if(type == "CHARGED"){
        return std::make_shared<CHARGEDparticleFilter>();
    } else if(type == "NEUTRAL"){
        return std::make_shared<NEUTRALparticleFilter>();
    } else if(type == "EM0"){
        return std::make_shared<EM0particleFilter>();
    } else if(type == "HAD0"){
        return std::make_shared<HAD0particleFilter>();
    } else if(type == "EM0HAD0"){
        return std::make_shared<EM0HAD0particleFilter>();
    } else if(type == "EM0CHARGED"){
        return std::make_shared<EM0CHARGEDparticleFilter>();
    } else if(type == "HAD0CHARGED"){
        return std::make_shared<HAD0CHARGEDparticleFilter>();
    } else {
        throw std::runtime_error("particleFilter::getParticleFilter: Unknown particleFilterType");
    }
}
