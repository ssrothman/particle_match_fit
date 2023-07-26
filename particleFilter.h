#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <vector>
#include <memory>
#include "SRothman/SimonTools/src/jets.h"

enum class particleFilterType{
    NONE = 0,
    ALL = 1,
    CHARGED = 2,
    NEUTRAL = 3,
    EM0 = 4,
    HAD0 = 5
};

class particleFilter{
public:
    particleFilter() {};
    virtual ~particleFilter(){};
    virtual bool pass(const particle& part) = 0;

    static std::shared_ptr<particleFilter> getParticleFilter(particleFilterType type);
};

#endif
