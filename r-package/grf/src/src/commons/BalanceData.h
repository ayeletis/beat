/*-------------------------------------------------------------------------------
  This file is part of generalized random forest (grf).

  grf is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  grf is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with grf. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef GRF_BALANCEDATA_H_
#define GRF_BALANCEDATA_H_

#include "commons/Data.h"
#include "Arma/rcpparma"
// [[Rcpp::depends(RcppArmadillo)]]

namespace grf
{

  /**
 * Data wrapper for GRF.
 * Serves as a read-only (immutable) wrapper of a column major (Fortran order)
 * array accessed through its pointer (data_ptr). This class does not own
 * data.
 *
 * The GRF data model is a contiguous array [X, Y, z, ...] of covariates X,
 * outcomes Y, and other optional variables z.
 *
 */
  class BalanceData : public Data
  {
  public:
    // set
    void set_target_avg_weights(arma::cube x);

    void set_target_weight_penalty(double target_weight_penalty);

    void set_target_weight_penalty_metric(std::string metric_type);

    void set_num_target_weight_cols(size_t num_cols);

    // get
    arma::vec get_target_weight_row(size_t x_var, size_t sample_id) const;

    double get_target_weight_penalty() const;

    std::string get_target_weight_penalty_metric() const;

    size_t get_num_target_weight_cols() const;

  protected:
    size_t num_target_weight_cols;

    double target_weight_penalty;
    std::string target_weight_penalty_metric;

    arma::cube target_avg_weights;
  };

} // namespace grf
#endif /* GRF_BALANCEDATA_H_ */
