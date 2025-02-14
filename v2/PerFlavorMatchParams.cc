#include "PerFlavorMatchParams.h"

matching::MatchParams::MatchParams(
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
        const std::string& flavor_filter_mode):
    dR_limiter(DeltaRLimiter::get_deltaRlimiter(
            dr_mode,
            dr_param1,
            dr_param2,
            dr_param3)),
    chi_sq_fn(
            ptres_mode,
            ptres_param1,
            ptres_param2,
            angres_mode,
            angres_param1,
            angres_param2,
            opp_charge_penalty,
            no_charge_penalty),
    charge_filter(ChargeFilter::get_charge_filter(charge_filter_mode)),
    flavor_filter(FlavorFilter::get_flavor_filter(flavor_filter_mode)) {}

matching::PerFlavorMatchParams::PerFlavorMatchParams() :
    ele_params(nullptr),
    mu_params(nullptr),
    hadch_params(nullptr),
    pho_params(nullptr),
    had0_params(nullptr) {}

void matching::PerFlavorMatchParams::setup_params(
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
        const std::string& flavor_filter_mode) {
    auto& target = get_target(flavor);
    if(target){
        throw std::invalid_argument("Flavor already set up");
    }

    target = std::make_unique<const MatchParams>(
            dr_mode,
            dr_param1,
            dr_param2,
            dr_param3,
            ptres_mode,
            ptres_param1,
            ptres_param2,
            angres_mode,
            angres_param1,
            angres_param2,
            opp_charge_penalty,
            no_charge_penalty,
            charge_filter_mode,
            flavor_filter_mode);
}

const matching::MatchParams& matching::PerFlavorMatchParams::get_params(Flavor flavor) const {
    auto& target = get_target(flavor);
    if(!target){
        printf("TRYING TO GET FLAVOR %d\n", flavor);
        throw std::invalid_argument("Flavor not set up");
    }
    return *target;
}

const matching::MatchParams& matching::PerFlavorMatchParams::get_params(const simon::particle& recopart) const {
    Flavor flavor;
    if(recopart.pdgid == 11){
        flavor = ELE;
    } else if(recopart.pdgid == 13){
        flavor = MU;
    } else if(recopart.pdgid == 22){
        flavor = PHO;
    } else if(recopart.charge !=0){
        flavor = HADCH;
    } else {
        flavor = HAD0;
    } 
    return get_params(flavor);
}

matching::MatchParamsPtr& matching::PerFlavorMatchParams::get_target(Flavor flavor) {
    switch(flavor){
        case ELE:
            return ele_params;
        case MU:
            return mu_params;
        case HADCH:
            return hadch_params;
        case PHO:
            return pho_params;
        case HAD0:
            return had0_params;
        default:
            throw std::invalid_argument("Invalid flavor");
    }
}

const matching::MatchParamsPtr& matching::PerFlavorMatchParams::get_target(Flavor flavor) const {
    switch(flavor){
        case ELE:
            return ele_params;
        case MU:
            return mu_params;
        case HADCH:
            return hadch_params;
        case PHO:
            return pho_params;
        case HAD0:
            return had0_params;
        default:
            throw std::invalid_argument("Invalid flavor");
    }
}

void matching::PerFlavorMatchParams::print_status() const {
    printf("PerFlavorMatchParams::status():\n");
    printf("\tELE: %s\n", ele_params ? "set" : "not set");
    printf("\tMU: %s\n", mu_params ? "set" : "not set");
    printf("\tHADCH: %s\n", hadch_params ? "set" : "not set");
    printf("\tPHO: %s\n", pho_params ? "set" : "not set");
    printf("\tHAD0: %s\n", had0_params ? "set" : "not set");
}


#ifdef CMSSW_GIT_HASH

matching::MatchParams::MatchParams(
        const edm::ParameterSet& params):
    MatchParams(
            params.getParameter<std::string>("dr_mode"),
            params.getParameter<double>("dr_param1"),
            params.getParameter<double>("dr_param2"),
            params.getParameter<double>("dr_param3"),
            params.getParameter<std::string>("ptres_mode"),
            params.getParameter<double>("ptres_param1"),
            params.getParameter<double>("ptres_param2"),
            params.getParameter<std::string>("angres_mode"),
            params.getParameter<double>("angres_param1"),
            params.getParameter<double>("angres_param2"),
            params.getParameter<double>("opp_charge_penalty"),
            params.getParameter<double>("no_charge_penalty"),
            params.getParameter<std::string>("charge_filter_mode"),
            params.getParameter<std::string>("flavor_filter_mode")) {}

void matching::MatchParams::fillPSetDescription(edm::ParameterSetDescription& desc) {
    desc.add<std::string>("dr_mode");
    desc.add<double>("dr_param1");
    desc.add<double>("dr_param2");
    desc.add<double>("dr_param3");
    desc.add<std::string>("ptres_mode");
    desc.add<double>("ptres_param1");
    desc.add<double>("ptres_param2");
    desc.add<std::string>("angres_mode");
    desc.add<double>("angres_param1");
    desc.add<double>("angres_param2");
    desc.add<double>("opp_charge_penalty");
    desc.add<double>("no_charge_penalty");
    desc.add<std::string>("charge_filter_mode");
    desc.add<std::string>("flavor_filter_mode");
}

void matching::PerFlavorMatchParams::setup_params(
        Flavor flavor,
        const edm::ParameterSet& params) {
    auto& target = get_target(flavor);
    if(target){
        throw std::invalid_argument("Flavor already set up");
    }
    target = std::make_unique<const MatchParams>(params);
}

#endif
