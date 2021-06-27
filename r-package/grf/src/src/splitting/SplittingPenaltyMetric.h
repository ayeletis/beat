/*
metrics for target weight splitting 
used in balanced instrumental and regression splitting 
*/
#include <RcppEigen.h>
#include "Arma/rcpparma"
// [[Rcpp::depends(RcppArmadillo)]]
namespace grf
{

    double calculate_target_weight_penalty(const double &target_weight_penalty_rate,
                                           const double &decrease_left,
                                           const double &decrease_right,
                                           const arma::vec &target_weights_avg_left,
                                           const arma::vec &target_weights_avg_right,
                                           const std::string &target_weight_penalty_metric);

}