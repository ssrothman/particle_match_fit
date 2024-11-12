#include "matchingUtil.h"

Eigen::MatrixXd fullmat(const unsigned nrow, const unsigned ncol,
                  const std::vector<std::pair<unsigned, unsigned>>& locs,
                  const std::vector<double>& vals){
    Eigen::MatrixXd ans(nrow, ncol);
    for(unsigned i = 0; i < locs.size(); i++){
        ans(locs[i].first, locs[i].second) = vals[i];
    }

    Eigen::VectorXd sums = ans.colwise().sum();
    for (unsigned i=0; i<sums.size(); ++i){
        if (sums[i] == 0){
            sums[i] = 1;
        }
    }
    ans.array().rowwise() /= sums.transpose().array();

    return ans;
}
