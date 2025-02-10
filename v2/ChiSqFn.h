#ifndef SROTHMAN_MATCHING_V2_CHISQFN_H
#define SROTHMAN_MATCHING_V2_CHISQFN_H

#include "ResFunc.h"
#include <string>

namespace matching {
    class ChiSqFn {
    public:
        ChiSqFn(const std::string& ptres_mode,
                const double ptres_param1,
                const double ptres_param2,
                const std::string& angres_mode,
                const double angres_param1,
                const double angres_param2,
                const double opp_charge_penalty,
                const double no_charge_penalty);

        double evaluate(const double pt1, 
                        const double eta1, 
                        const double phi1, 
                        const int charge1,
                        const double pt2, 
                        const double eta2, 
                        const double phi2,
                        const int charge2) const;

    private:
        const ResFuncPtr ptresfunc, angresfunc;
        const double opp_charge_penalty, no_charge_penalty;
    };
};


#endif
