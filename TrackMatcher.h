#ifndef SROTHMAN_MATCHING_V2_TRACKMATCHER_H
#define SROTHMAN_MATCHING_V2_TRACKMATCHER_H

#include "SRothman/SimonTools/src/jet.h"
#include "PerFlavorMatchParams.h"

#include <string>
#include <vector>

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace matching {
    struct matchidxs {
        size_t iReco, iGen;
        matchidxs(size_t iReco, size_t iGen) : iReco(iReco), iGen(iGen) {}
    };
    using matchvec = std::vector<matchidxs>;

    class TrackMatcher {
    public:
        TrackMatcher(
                //jet parameters
                const double jet_dR_threshold,

                //global params
                const double max_chisq,

                //electron params
                const std::string& ele_dr_mode,
                const double ele_dr_param1,
                const double ele_dr_param2,
                const double ele_dr_param3,
                const std::string& ele_ptres_mode,
                const double ele_ptres_param1,
                const double ele_ptres_param2,
                const std::string& ele_angres_mode,
                const double ele_angres_param1,
                const double ele_angres_param2,
                const double ele_opp_charge_penalty,
                const double ele_no_charge_penalty,
                const std::string& ele_charge_filter_mode,
                const std::string& ele_flavor_filter_mode,

                //muon params
                const std::string& mu_dr_mode,
                const double mu_dr_param1,
                const double mu_dr_param2,
                const double mu_dr_param3,
                const std::string& mu_ptres_mode,
                const double mu_ptres_param1,
                const double mu_ptres_param2,
                const std::string& mu_angres_mode,
                const double mu_angres_param1,
                const double mu_angres_param2,
                const double mu_opp_charge_penalty,
                const double mu_no_charge_penalty,
                const std::string& mu_charge_filter_mode,
                const std::string& mu_flavor_filter_mode,

                //charged hadron params
                const std::string& hadch_dr_mode,
                const double hadch_dr_param1,
                const double hadch_dr_param2,
                const double hadch_dr_param3,
                const std::string& hadch_ptres_mode,
                const double hadch_ptres_param1,
                const double hadch_ptres_param2,
                const std::string& hadch_angres_mode,
                const double hadch_angres_param1,
                const double hadch_angres_param2,
                const double hadch_opp_charge_penalty,
                const double hadch_no_charge_penalty,
                const std::string& hadch_charge_filter_mode,
                const std::string& hadch_flavor_filter_mode);

        void matchJets(
            const std::vector<simon::jet>& recojets,
            const std::vector<simon::jet>& genjets,
            matchvec& matches);

        void matchParticles(
            const simon::jet& recojet,
            const simon::jet& genjet,
            Eigen::MatrixXd& tmat);

#ifdef CMSSW_GIT_HASH
        TrackMatcher(const edm::ParameterSet& iConfig);

        static void fillPSetDescription(edm::ParameterSetDescription& desc);
#endif

    private:
        const double jet_dR_threshold;
        const double max_chisq;

        PerFlavorMatchParams particle_params;
    };
};

#endif
