/*
metrics for target weight splitting 
used in balanced instrumental and regression splitting 
*/
#include <RcppEigen.h>
 #include "Arma/rcpparma" 
// [[Rcpp::depends(RcppArmadillo)]]
namespace grf
{

    double calculate_target_weight_penalty(double target_weight_penalty_rate,
                                           double decrease_left,
                                           double decrease_right,
                                           arma::rowvec target_weights_avg_left,
                                           arma::rowvec target_weights_avg_right,
                                           std::string target_weight_penalty_metric);

}