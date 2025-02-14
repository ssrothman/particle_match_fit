#include "TrackMatcher.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include <numeric>

static constexpr double INF = std::numeric_limits<double>::infinity();

matching::TrackMatcher::TrackMatcher(
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
        const std::string& hadch_flavor_filter_mode):

    jet_dR_threshold(jet_dR_threshold),
    max_chisq(max_chisq),
    particle_params() {
    
    particle_params.setup_params(
        PerFlavorMatchParams::ELE,
        ele_dr_mode,
        ele_dr_param1,
        ele_dr_param2,
        ele_dr_param3,
        ele_ptres_mode,
        ele_ptres_param1,
        ele_ptres_param2,
        ele_angres_mode,
        ele_angres_param1,
        ele_angres_param2,
        ele_opp_charge_penalty,
        ele_no_charge_penalty,
        ele_charge_filter_mode,
        ele_flavor_filter_mode);

    particle_params.setup_params(
        PerFlavorMatchParams::MU,
        mu_dr_mode,
        mu_dr_param1,
        mu_dr_param2,
        mu_dr_param3,
        mu_ptres_mode,
        mu_ptres_param1,
        mu_ptres_param2,
        mu_angres_mode,
        mu_angres_param1,
        mu_angres_param2,
        mu_opp_charge_penalty,
        mu_no_charge_penalty,
        mu_charge_filter_mode,
        mu_flavor_filter_mode);

    particle_params.setup_params(
        PerFlavorMatchParams::HADCH,
        hadch_dr_mode,
        hadch_dr_param1,
        hadch_dr_param2,
        hadch_dr_param3,
        hadch_ptres_mode,
        hadch_ptres_param1,
        hadch_ptres_param2,
        hadch_angres_mode,
        hadch_angres_param1,
        hadch_angres_param2,
        hadch_opp_charge_penalty,
        hadch_no_charge_penalty,
        hadch_charge_filter_mode,
        hadch_flavor_filter_mode);
}

template <typename T>
static void match_one_to_one(
        const std::vector<T>& recovec,
        const std::vector<T>& genvec,
        const matching::PerFlavorMatchParams& particle_params,
        const double max_chisq,
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
        const auto& theparms = particle_params.get_params(reco);
        
        double best_chisq = INF;
        int best_igen = -1;

        for(size_t iGen : gen_ptorder){
            if(gen_used[iGen]) continue;

            const auto& gen = genvec[iGen];

            double dR = simon::deltaR(gen.eta, gen.phi,
                                     reco.eta, reco.phi);

            double dRlim = theparms.dR_limiter->evaluate(
                    reco.pt, reco.eta, reco.phi);
            if(dR > dRlim) continue;

            if(!theparms.charge_filter->evaluate(
                    reco.charge, gen.charge)) continue;

            if(!theparms.flavor_filter->evaluate(
                    gen.charge,
                    gen.pdgid)) continue;

            double chisq = theparms.chi_sq_fn.evaluate(
                    reco.pt, reco.eta, 
                    reco.phi, reco.charge,
                    gen.pt, gen.eta, 
                    gen.phi, gen.charge);

            if(chisq < best_chisq){
                best_chisq = chisq;
                best_igen = iGen;
            }
        }//end gen loop
        if(best_igen>=0 && best_chisq < max_chisq){
            matches.emplace_back(iReco, best_igen);
        }
    }//end reco loop
}//end match_one_to_one()

void matching::TrackMatcher::matchJets(
        const std::vector<simon::jet>& recojets,
        const std::vector<simon::jet>& genjets,
        matchvec& matches){
    matches.clear();

    std::vector<size_t> gen_ptorder(genjets.size());
    std::iota(gen_ptorder.begin(), gen_ptorder.end(), 0);
    std::sort(gen_ptorder.begin(), gen_ptorder.end(),
            [&](size_t i1, size_t i2){
                return genjets[i1].pt > genjets[i2].pt;
            });

    std::vector<size_t> reco_ptorder(recojets.size());
    std::iota(reco_ptorder.begin(), reco_ptorder.end(), 0);
    std::sort(reco_ptorder.begin(), reco_ptorder.end(),
            [&](size_t i1, size_t i2){
                return recojets[i1].pt > recojets[i2].pt;
            });

    std::vector<bool> gen_used(genjets.size(), false);
    for(const size_t iRecoJet : reco_ptorder){
        double best_dR = INF;
        int matched_gen = -1;
        const auto& recojet = recojets[iRecoJet];
        for(const size_t iGenJet : gen_ptorder){
            if(gen_used[iGenJet]) continue;

            const auto& genjet = genjets[iGenJet];

            double dR = simon::deltaR(genjet.eta, genjet.phi,
                                     recojet.eta, recojet.phi);

            if(dR < best_dR){
                best_dR = dR;
                matched_gen = iGenJet;
            }
        }//end gen loop
        if(best_dR < jet_dR_threshold){
            matches.emplace_back(iRecoJet, matched_gen);
            gen_used[matched_gen] = true;
        }
    }
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
    match_one_to_one(
            recoparts, genparts,
            particle_params,
            max_chisq,
            matches);

    for(const auto& match : matches){
        tmat(match.iReco, match.iGen) = 1;
    }
}

#ifdef CMSSW_GIT_HASH
matching::TrackMatcher::TrackMatcher(const edm::ParameterSet& iConfig) :
    jet_dR_threshold(iConfig.getParameter<double>("jet_dR_threshold")),
    max_chisq(iConfig.getParameter<double>("max_chisq")),
    particle_params() {

    particle_params.setup_params(
        PerFlavorMatchParams::ELE,
        iConfig.getParameter<edm::ParameterSet>("Electrons")
    );
    particle_params.setup_params(
        PerFlavorMatchParams::MU,
        iConfig.getParameter<edm::ParameterSet>("Muons")
    );
    particle_params.setup_params(
        PerFlavorMatchParams::HADCH,
        iConfig.getParameter<edm::ParameterSet>("ChargedHadrons")
    );
}

void matching::TrackMatcher::fillPSetDescription(edm::ParameterSetDescription& desc){
    desc.add<double>("jet_dR_threshold");
    desc.add<double>("max_chisq");

    edm::ParameterSetDescription ele_desc;
    MatchParams::fillPSetDescription(ele_desc);
    desc.add<edm::ParameterSetDescription>("Electrons", ele_desc);

    edm::ParameterSetDescription mu_desc;
    MatchParams::fillPSetDescription(mu_desc);
    desc.add<edm::ParameterSetDescription>("Muons", mu_desc);

    edm::ParameterSetDescription hadch_desc;
    MatchParams::fillPSetDescription(hadch_desc);
    desc.add<edm::ParameterSetDescription>("ChargedHadrons", hadch_desc);
}
#endif
