#ifndef REFINE_PREFIT_H
#define REFINE_PREFIT_H

#include <vector>
#include <memory>
#include "SRothman/SimonTools/src/jet.h"
#include <unordered_map>
#include <string>

class prefitRefiner{
    public:
        prefitRefiner() {}
        virtual ~prefitRefiner() {}

        typedef std::unordered_map<unsigned, std::vector<unsigned>> matchMap;

        virtual matchMap refine(const matchMap& recoToGen,
                                const simon::jet& recojet,
                                const simon::jet& genjet) const = 0;
                
        static std::shared_ptr<prefitRefiner> get(const std::string& behavior);

    protected:
        matchMap invertMap(const matchMap& m) const;

};
#endif
