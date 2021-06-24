/*
metrics for target weight splitting 
used in balanced instrumental and regression splitting 
*/
#include <RcppEigen.h>

namespace grf
{

    double calculate_target_weight_penalty(double target_weight_penalty_rate,
                                           double decrease_left,
                                           double decrease_right,
                                           Eigen::VectorXd target_weights_avg_left,
                                           Eigen::VectorXd target_weights_avg_right,
                                           std::string target_weight_penalty_metric);

}