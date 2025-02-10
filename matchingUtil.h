#ifndef MATCHING_UTIL_H
#define MATCHING_UTIL_H

#include <vector>
#include "SRothman/SimonTools/src/jet.h"
#include "SRothman/SimonTools/src/deltaR.h"
#include "SRothman/SimonTools/src/util.h"
#include <string>

template<typename T>
T chisquared(const T& recopt, const T& recoeta, const T& recophi,
             const T& genpt, const T& geneta, const T& genphi,
             const T& dpt, const T& deta, const T& dphi,
             bool wpt){
    double errpt = (recopt - genpt) / dpt;
    double erreta = (recoeta - geneta) / deta;
    double errphi = simon::deltaPhi(recophi, genphi) / dphi;

    if (wpt){
        return errpt*errpt + erreta*erreta + errphi*errphi;
    } else {
        return erreta*erreta + errphi*errphi;
    }
}

inline double chisquared(const simon::particle& reco, const simon::particle& gen,
                         bool wpt){
    return chisquared(reco.pt, reco.eta, reco.phi, 
                      gen.pt, gen.eta, gen.phi,
                      reco.dpt, reco.deta, reco.dphi,
                      wpt);
}

Eigen::MatrixXd fullmat(const unsigned nrow, const unsigned ncol,
                  const std::vector<std::pair<unsigned, unsigned>>& locs,
                  const std::vector<double>& vals);

#endif
