/*
metrics for target weight splitting
used in balanced instrumental and regression splitting
*/
#include <RcppEigen.h>
#include "SplittingPenaltyMetric.h"

namespace grf
{

    double calculate_target_weight_penalty(double target_weight_penalty_rate,
                                           double decrease_left,
                                           double decrease_right,
                                           Eigen::VectorXd target_weights_avg_left,
                                           Eigen::VectorXd target_weights_avg_right,
                                           std::string target_weight_penalty_metric)
    // return penalty for target weight imbalance
    {
        double imbalance = 0;

        // rate methods
        if (target_weight_penalty_metric == "split_l2_norm_rate")
        {
            // std::cout << "var" << var << "decrease:" << decrease << "penalty:" << penalty_target_weight << "\n";
            imbalance = target_weights_avg_left.lpNorm<2>() * decrease_left + target_weights_avg_right.lpNorm<2>() * decrease_right;
        }
        else if (target_weight_penalty_metric == "euclidean_distance_rate")
        {
            imbalance = sqrt((target_weights_avg_left - target_weights_avg_right).pow(2).sum()); //>=0
            imbalance = (decrease_left + decrease_right) * imbalance;
        }
        else if (target_weight_penalty_metric == "cosine_similarity_rate")
        {
            double upper = (target_weights_avg_left * target_weights_avg_right).sum();
            double lower = sqrt(target_weights_avg_left.pow(2).sum()) * sqrt(target_weights_avg_right.pow(2).sum());
            double cosine_similarity = upper / lower; //−1 = exactly opposite, 1 = exactly the same
            imbalance = 1 - cosine_similarity;        // in [0, 2]
            imbalance = (decrease_left + decrease_right) * imbalance;
        }

        // direct methods
        else if (target_weight_penalty_metric == "split_l2_norm")
        {
            // std::cout << "var" << var << "decrease:" << decrease << "penalty:" << penalty_target_weight << "\n";
            imbalance = target_weights_avg_left.lpNorm<2>() + target_weights_avg_right.lpNorm<2>();
        }
        else if (target_weight_penalty_metric == "euclidean_distance")
        {
            imbalance = sqrt((target_weights_avg_left - target_weights_avg_right).pow(2).sum()); //>=0
        }
        else if (target_weight_penalty_metric == "cosine_similarity")
        {
            double upper = (target_weights_avg_left * target_weights_avg_right).sum();
            double lower = sqrt(target_weights_avg_left.pow(2).sum()) * sqrt(target_weights_avg_right.pow(2).sum());
            double cosine_similarity = upper / lower; //−1 = exactly opposite, 1 = exactly the same
            imbalance = 1 - cosine_similarity;        // in [0, 2]
        }
        // ends
        else
        { // the R wrapper function shall have filtered out
            throw std::invalid_argument("penalty metric type is not available");
        }

        return imbalance * target_weight_penalty_rate;
    }

} // namespace grf