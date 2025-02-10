#ifndef SROTHMAN_MATCHING_V2_TRACKMATCHER_H
#define SROTHMAN_MATCHING_V2_TRACKMATCHER_H

#include "SRothman/SimonTools/src/jet.h"
#include "ChiSqFn.h"
#include "DeltaRLimiter.h"

#include <string>
#include <vector>

namespace matching {
    class TrackMatcher {
    public:
        TrackMatcher(const std::string& jet_dr_mode,
                     const double jet_dr_param1,
                     const double jet_dr_param2,
                     const double jet_dr_param3,
                     const std::string& jet_ptres_mode,
                     const double jet_ptres_param1,
                     const double jet_ptres_param2,
                     const std::string& jet_angres_mode,
                     const double jet_angres_param1,
                     const double jet_angres_param2,
                     const std::string& particle_dr_mode,
                     const double particle_dr_param1,
                     const double particle_dr_param2,
                     const double particle_dr_param3,
                     const std::string& particle_ptres_mode,
                     const double particle_ptres_param1,
                     const double particle_ptres_param2,
                     const std::string& particle_angres_mode,
                     const double particle_angres_param1,
                     const double particle_angres_param2,
                     const double opp_charge_penalty,
                     const double no_charge_penalty);

        void matchJets(
            const std::vector<simon::jet>& genjets,
            const std::vector<simon::jet>& recojets,
            std::vector<std::pair<simon::jet const*, simon::jet const *>>& matches);

        void matchParticles(
            const simon::jet& genjet,
            const simon::jet& recojet,
            Eigen::MatrixXd& tmat);

    private:
        const DeltaRLimiterPtr jet_dR_limiter;
        const ChiSqFn jet_chisq_fn;

        const DeltaRLimiterPtr particle_dR_limiter;
        const ChiSqFn particle_chisq_fn;
    };
};

#endif
