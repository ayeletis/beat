/*
metrics for target weight splitting
used in balanced instrumental and regression splitting
*/
#include <RcppEigen.h>
#include "SplittingPenaltyMetric.h"
#include "Arma/rcpparma"
// [[Rcpp::depends(RcppArmadillo)]]

namespace grf
{

    double calculate_target_weight_penalty(const double &target_weight_penalty_rate,
                                           const double &decrease_left,
                                           const double &decrease_right,
                                           const arma::vec &target_weight_avg_left,
                                           const arma::vec &target_weight_avg_right,
                                           const std::string &target_weight_penalty_metric)
    // return penalty for target weight imbalance
    {
        double imbalance = 0;
        // arma::vec target_left = arma::vectorise(target_weight_avg_left);
        // arma::vec target_right = arma::vectorise(target_weight_avg_right);

        // rate methods
        if (target_weight_penalty_metric == "split_l2_norm_rate")
        {
            //            double left = sqrt(arma::sum(arma::square(target_left)));
            //            double right = sqrt(arma::sum(arma::square(target_right)));
            //            imbalance = left * decrease_left + right * decrease_right;
            imbalance = arma::norm(target_weight_avg_left, 2) * decrease_left + arma::norm(target_weight_avg_right, 2) * decrease_right;
        }
        else if (target_weight_penalty_metric == "euclidean_distance_rate")
        {

            imbalance = arma::norm((target_weight_avg_left - target_weight_avg_right), 2);
            imbalance = (decrease_left + decrease_right) * imbalance;
        }
        else if (target_weight_penalty_metric == "cosine_similarity_rate")
        {

            double upper = arma::sum((target_weight_avg_left % target_weight_avg_right));
            double lower = arma::norm(target_weight_avg_left, 2) * arma::norm(target_weight_avg_right, 2);

            double cosine_similarity = upper / lower; //−1 = exactly opposite, 1 = exactly the same
            imbalance = 1 - cosine_similarity / 2;    // in [0, 1]
            imbalance = (decrease_left + decrease_right) * imbalance;
        }

        // direct methods
        else if (target_weight_penalty_metric == "split_l2_norm")
        {
            imbalance = arma::norm(target_weight_avg_left, 2) + arma::norm(target_weight_avg_right, 2);
        }
        else if (target_weight_penalty_metric == "euclidean_distance")
        {
            imbalance = arma::norm((target_weight_avg_left - target_weight_avg_right), 2);
        }
        else if (target_weight_penalty_metric == "cosine_similarity")
        {

            double upper = arma::sum((target_weight_avg_left % target_weight_avg_right));
            double lower = arma::norm(target_weight_avg_left, 2) * arma::norm(target_weight_avg_right, 2);

            double cosine_similarity = upper / lower; //−1 = exactly opposite, 1 = exactly the same
            imbalance = 1 - cosine_similarity / 2;    // in [0, 1]
        }
        // ends
        else
        { // the R wrapper function shall have filtered out
            throw std::invalid_argument("penalty metric type is not available");
        }

        return imbalance * target_weight_penalty_rate;
    }

} // namespace grf
