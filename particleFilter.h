#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <vector>
#include <memory>
#include "SRothman/SimonTools/src/jets.h"
#include <string>

class particleFilter{
public:
    particleFilter() {};
    virtual ~particleFilter(){};
    virtual bool pass(const particle& part) = 0;

    static std::shared_ptr<particleFilter> getParticleFilter(const std::string& type);
};

#endif
