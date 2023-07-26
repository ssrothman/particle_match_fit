#include "matchingUtil.h"

arma::mat fullmat(const unsigned nrow, const unsigned ncol,
                  const std::vector<std::pair<unsigned, unsigned>>& locs,
                  const std::vector<double>& vals){
    arma::mat ans(nrow, ncol, arma::fill::zeros);
    for(unsigned i = 0; i < locs.size(); i++){
        ans(locs[i].first, locs[i].second) = vals[i];
    }

    arma::rowvec sums = arma::sum(ans, 0);
    sums.replace(0, 1);
    ans.each_row() /= sums;

    return ans;
}
