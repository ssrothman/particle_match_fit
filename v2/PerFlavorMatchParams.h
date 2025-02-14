#ifndef SROTHMAN_MATCHING_PERFLAVORMATCHPARAMS_H
#define SROTHMAN_MATCHING_PERFLAVORMATCHPARAMS_H

#include "DeltaRLimiter.h"
#include "ChiSqFn.h"
#include "ChargeFilter.h"
#include "FlavorFilter.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace matching{
    class MatchParams {
    public:
        MatchParams(
            //dR_limiter params
            const std::string& dr_mode,
            const double dr_param1,
            const double dr_param2,
            const double dr_param3,
            //chi_sq_fn params
            const std::string& ptres_mode,
            const double ptres_param1,
            const double ptres_param2,
            const std::string& angres_mode,
            const double angres_param1,
            const double angres_param2,
            const double opp_charge_penalty,
            const double no_charge_penalty,
            //charge_filter params
            const std::string& charge_filter_mode,
            //flavor_filter params
            const std::string& flavor_filter_mode);

        const DeltaRLimiterPtr dR_limiter;
        const ChiSqFn chi_sq_fn;
        const ChargeFilterPtr charge_filter;
        const FlavorFilterPtr flavor_filter;

#ifdef CMSSW_GIT_HASH
        MatchParams(
            const edm::ParameterSet& params);

        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif
    };

    using MatchParamsPtr = std::unique_ptr<const MatchParams>;

    class PerFlavorMatchParams {
    public:
        PerFlavorMatchParams();

        enum Flavor{
            ELE=0,
            MU=1,
            HADCH=2,
            PHO=3,
            HAD0=4
        };
        void setup_params(
                Flavor flavor,
                //dR_limiter params
                const std::string& dr_mode,
                const double dr_param1,
                const double dr_param2,
                const double dr_param3,
                //chi_sq_fn params
                const std::string& ptres_mode,
                const double ptres_param1,
                const double ptres_param2,
                const std::string& angres_mode,
                const double angres_param1,
                const double angres_param2,
                const double opp_charge_penalty,
                const double no_charge_penalty,
                //charge_filter params
                const std::string& charge_filter_mode,
                //flavor_filter params
                const std::string& flavor_filter_mode);
    
#ifdef CMSSW_GIT_HASH
        void setup_params(
                Flavor flavor,
                const edm::ParameterSet& params);
#endif

        const MatchParams& get_params(Flavor flavor) const;
        const MatchParams& get_params(const simon::particle& recopart) const;

        void print_status() const;

    private:
        MatchParamsPtr ele_params;
        MatchParamsPtr mu_params;
        MatchParamsPtr hadch_params;
        MatchParamsPtr pho_params;
        MatchParamsPtr had0_params;

        MatchParamsPtr& get_target(Flavor flavor);
        const MatchParamsPtr& get_target(Flavor flavor) const;
    };
};

#endif
