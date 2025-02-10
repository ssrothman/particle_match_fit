#include "TrackMatcher.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include <numeric>

static constexpr double INF = std::numeric_limits<double>::infinity();

matching::TrackMatcher::TrackMatcher(
        const std::string& jet_dr_mode,
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
        const double no_charge_penalty):
    jet_dR_limiter(DeltaRLimiter::get_deltaRlimiter(
            jet_dr_mode,
            jet_dr_param1,
            jet_dr_param2,
            jet_dr_param3)),
    jet_chisq_fn(
            jet_ptres_mode,
            jet_ptres_param1,
            jet_ptres_param2,
            jet_angres_mode,
            jet_angres_param1,
            jet_angres_param2,
            0.0, 0.0),
    particle_dR_limiter(DeltaRLimiter::get_deltaRlimiter(
            particle_dr_mode,
            particle_dr_param1,
            particle_dr_param2,
            particle_dr_param3)),
    particle_chisq_fn(particle_ptres_mode,
            particle_ptres_param1,
            particle_ptres_param2,
            particle_angres_mode,
            particle_angres_param1,
            particle_angres_param2,
            opp_charge_penalty,
            no_charge_penalty) {}

template <typename T>
static void match_one_to_one(
        const std::vector<T>& genvec,
        const std::vector<T>& recovec,
        const matching::DeltaRLimiterPtr& dRlimiter,
        const matching::ChiSqFn& chisq_fn,
        std::vector<std::pair<size_t, size_t>>& matches){

    matches.clear();

    std::vector<size_t> gen_ptorder(genvec.size());
    std::iota(gen_ptorder.begin(), gen_ptorder.end(), 0);
    std::sort(gen_ptorder.begin(), gen_ptorder.end(),
            [&](size_t i1, size_t i2){
                return genvec[i1].pt > genvec[i2].pt;
            });

    std::vector<size_t> reco_ptorder(recovec.size());
    std::iota(reco_ptorder.begin(), reco_ptorder.end(), 0);
    std::sort(reco_ptorder.begin(), reco_ptorder.end(),
            [&](size_t i1, size_t i2){
                return recovec[i1].pt > recovec[i2].pt;
            });

    std::vector<bool> gen_used(genvec.size(), false);

    for(size_t iReco : reco_ptorder){
        const auto& reco = recovec[iReco];
        printf("reco pt: %f\n", reco.pt);
        
        double best_chisq = INF;
        int best_igen = -1;

        for(size_t iGen : gen_ptorder){
            if(gen_used[iGen]) continue;

            const auto& gen = genvec[iGen];

            double dR = simon::deltaR(gen.eta, gen.phi,
                                     reco.eta, reco.phi);
            double dRlim = dRlimiter->evaluate(
                    reco.pt, reco.eta, reco.phi);
            if(dR > dRlim) continue;

            double chisq = chisq_fn.evaluate(
                    reco.pt, reco.eta, 
                    reco.phi, reco.charge,
                    gen.pt, gen.eta, 
                    gen.phi, gen.charge);

            if(chisq < best_chisq){
                best_chisq = chisq;
                best_igen = iGen;
            }
        }//end gen loop
        matches.emplace_back(iReco, best_igen);
    }//end reco loop
}//end match_one_to_one()

void matching::TrackMatcher::matchJets(
        const std::vector<simon::jet>& genjets,
        const std::vector<simon::jet>& recojets,
        std::vector<std::pair<simon::jet const *, simon::jet const *>>& matches){
    matches.clear();

    std::vector<std::pair<size_t, size_t>> match_indices;

    match_one_to_one(genjets, recojets,
                     jet_dR_limiter, jet_chisq_fn, 
                     match_indices);

    for(const auto& match : match_indices){
        matches.emplace_back(
            const_cast<simon::jet*>(&genjets[match.second]),
            const_cast<simon::jet*>(&recojets[match.first]));
    }
}

void matching::TrackMatcher::matchParticles(
        const simon::jet& genjet,
        const simon::jet& recojet,
        Eigen::MatrixXd& tmat){

    tmat.resize(recojet.nPart, genjet.nPart);

    const auto& genparts = genjet.particles;
    const auto& recoparts = recojet.particles;

    std::vector<std::pair<size_t, size_t>> matches;
    match_one_to_one(genparts, recoparts,
                     particle_dR_limiter, particle_chisq_fn,
                     matches);

    for(const auto& match : matches){
        tmat(match.first, match.second) = 1;
    }
}
