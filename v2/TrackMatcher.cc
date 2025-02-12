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
        const std::vector<T>& recovec,
        const std::vector<T>& genvec,
        const matching::DeltaRLimiterPtr& dRlimiter,
        const matching::ChiSqFn& chisq_fn,
        matching::matchvec& matches){

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
        if(best_igen>=0){
            matches.emplace_back(iReco, best_igen);
        }
    }//end reco loop
}//end match_one_to_one()

void matching::TrackMatcher::matchJets(
        const std::vector<simon::jet>& recojets,
        const std::vector<simon::jet>& genjets,
        matchvec& matches){
    matches.clear();

    match_one_to_one(recojets, genjets,
                     jet_dR_limiter, jet_chisq_fn, 
                     matches);
}

void matching::TrackMatcher::matchParticles(
        const simon::jet& recojet,
        const simon::jet& genjet,
        Eigen::MatrixXd& tmat){

    tmat.resize(recojet.nPart, genjet.nPart);
    tmat.setZero();

    const auto& genparts = genjet.particles;
    const auto& recoparts = recojet.particles;

    matchvec matches;
    match_one_to_one(recoparts, genparts,
                     particle_dR_limiter, particle_chisq_fn,
                     matches);

    for(const auto& match : matches){
        tmat(match.iReco, match.iGen) = 1;
    }
}

#ifdef CMSSW_GIT_HASH
matching::TrackMatcher::TrackMatcher(const edm::ParameterSet& iConfig) :
    TrackMatcher(
            iConfig.getParameter<std::string>("jet_dr_mode"),
            iConfig.getParameter<double>("jet_dr_param1"),
            iConfig.getParameter<double>("jet_dr_param2"),
            iConfig.getParameter<double>("jet_dr_param3"),
            iConfig.getParameter<std::string>("jet_ptres_mode"),
            iConfig.getParameter<double>("jet_ptres_param1"),
            iConfig.getParameter<double>("jet_ptres_param2"),
            iConfig.getParameter<std::string>("jet_angres_mode"),
            iConfig.getParameter<double>("jet_angres_param1"),
            iConfig.getParameter<double>("jet_angres_param2"),
            iConfig.getParameter<std::string>("particle_dr_mode"),
            iConfig.getParameter<double>("particle_dr_param1"),
            iConfig.getParameter<double>("particle_dr_param2"),
            iConfig.getParameter<double>("particle_dr_param3"),
            iConfig.getParameter<std::string>("particle_ptres_mode"),
            iConfig.getParameter<double>("particle_ptres_param1"),
            iConfig.getParameter<double>("particle_ptres_param2"),
            iConfig.getParameter<std::string>("particle_angres_mode"),
            iConfig.getParameter<double>("particle_angres_param1"),
            iConfig.getParameter<double>("particle_angres_param2"),
            iConfig.getParameter<double>("opp_charge_penalty"),
            iConfig.getParameter<double>("no_charge_penalty")) {}

void matching::TrackMatcher::fillPSetDescription(edm::ParameterSetDescription& desc){
    desc.add<std::string>("jet_dr_mode");
    desc.add<double>("jet_dr_param1");
    desc.add<double>("jet_dr_param2");
    desc.add<double>("jet_dr_param3");
    desc.add<std::string>("jet_ptres_mode");
    desc.add<double>("jet_ptres_param1");
    desc.add<double>("jet_ptres_param2");
    desc.add<std::string>("jet_angres_mode");
    desc.add<double>("jet_angres_param1");
    desc.add<double>("jet_angres_param2");
    desc.add<std::string>("particle_dr_mode");
    desc.add<double>("particle_dr_param1");
    desc.add<double>("particle_dr_param2");
    desc.add<double>("particle_dr_param3");
    desc.add<std::string>("particle_ptres_mode");
    desc.add<double>("particle_ptres_param1");
    desc.add<double>("particle_ptres_param2");
    desc.add<std::string>("particle_angres_mode");
    desc.add<double>("particle_angres_param1");
    desc.add<double>("particle_angres_param2");
    desc.add<double>("opp_charge_penalty");
    desc.add<double>("no_charge_penalty");
}
#endif
