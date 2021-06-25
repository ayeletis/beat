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

    double calculate_target_weight_penalty(double target_weight_penalty_rate,
                                           double decrease_left,
                                           double decrease_right,
                                           arma::rowvec target_weight_avg_left,
                                           arma::rowvec target_weight_avg_right,
                                           std::string target_weight_penalty_metric)
    // return penalty for target weight imbalance
    {
        double imbalance = 0;

        // rate methods
        if (target_weight_penalty_metric == "split_l2_norm_rate")
        {

            imbalance = arma::norm(target_weight_avg_left, 2) * decrease_left + arma::norm(target_weight_avg_right, 2) * decrease_right;
        }
        else if (target_weight_penalty_metric == "euclidean_distance_rate")
        {
            arma::rowvec gap = target_weight_avg_left - target_weight_avg_right;
            imbalance = sqrt(arma::sum(gap % gap));
            imbalance = (decrease_left + decrease_right) * imbalance;
        }
        else if (target_weight_penalty_metric == "cosine_similarity_rate")
        {
            double upper = arma::sum((target_weight_avg_left % target_weight_avg_right));
            double lower = sqrt(arma::sum(target_weight_avg_left % target_weight_avg_left)) * sqrt(arma::sum(target_weight_avg_right % target_weight_avg_right));

            double cosine_similarity = upper / lower; //−1 = exactly opposite, 1 = exactly the same
            imbalance = 1 - cosine_similarity / 2;    // in [0, 1]
            imbalance = (decrease_left + decrease_right) * imbalance;
        }

        // direct methods
        else if (target_weight_penalty_metric == "split_l2_norm")
        {
            // std::cout << "var" << var << "decrease:" << decrease << "penalty:" << penalty_target_weight << "\n";
            imbalance = arma::norm(target_weight_avg_left, 2) + arma::norm(target_weight_avg_right, 2);
        }
        else if (target_weight_penalty_metric == "euclidean_distance")

        {
            arma::rowvec gap = target_weight_avg_left - target_weight_avg_right;
            imbalance = sqrt(arma::sum(gap % gap));
        }
        else if (target_weight_penalty_metric == "cosine_similarity")
        {
            double upper = arma::sum((target_weight_avg_left % target_weight_avg_right));
            double lower = sqrt(arma::sum(target_weight_avg_left % target_weight_avg_left)) * sqrt(arma::sum(target_weight_avg_right % target_weight_avg_right));

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