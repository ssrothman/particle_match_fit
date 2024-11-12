#ifndef REFINE_PREFIT_H
#define REFINE_PREFIT_H

#include <vector>
#include <memory>
#include "SRothman/armadillo-12.2.0/include/armadillo"
#include "SRothman/SimonTools/src/jets.h"
#include <unordered_map>
#include <string>

class prefitRefiner{
    public:
        prefitRefiner() {}
        virtual ~prefitRefiner() {}

        typedef std::unordered_map<unsigned, std::vector<unsigned>> matchMap;

        virtual matchMap refine(const matchMap& recoToGen,
                                const jet& recojet,
                                const jet& genjet) const = 0;
                
        static std::shared_ptr<prefitRefiner> get(const std::string& behavior);

    protected:
        matchMap invertMap(const matchMap& m) const;

};
#endif
