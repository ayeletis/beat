# Explanation on the Source Code of Casual Forest Tree Splitting Algorithm 

A general description on the whole process can be found at [here](https://grf-labs.github.io/grf/REFERENCE.html#training). 

The R code of ``casual_forest`` is pointed to the c++ function  ``causal_train``, which is defined below. Note, the instrument variable is the same as the treatment variable (please see the parameter ``treatment_index`` ), and they are binary variables. The default value for sample weight is 1. 
```c++
//grf/src/CaualForestBindings.cpp
Rcpp::List causal_train(Rcpp::NumericMatrix train_matrix,
                        Eigen::SparseMatrix<double> sparse_train_matrix,
                        size_t outcome_index,
                        size_t treatment_index,
                        size_t sample_weight_index,
                        Rcpp::List target_avg_weights,
                        double target_weight_penalty,
                        bool use_sample_weights,
                        unsigned int mtry,
                        unsigned int num_trees,
                        unsigned int min_node_size,
                        double sample_fraction,
                        bool honesty,
                        double honesty_fraction,
                        bool honesty_prune_leaves,
                        size_t ci_group_size,
                        double reduced_form_weight,
                        double alpha,
                        double imbalance_penalty,
                        bool stabilize_splits,
                        std::vector<size_t> clusters,
                        unsigned int samples_per_cluster,
                        bool compute_oob_predictions,
                        unsigned int num_threads,
                        unsigned int seed)
{
    // set up relabelling, splitting, and prediction strategies
  ForestTrainer trainer = instrumental_trainer(reduced_form_weight, stabilize_splits);
 
 // set up the input data
  std::unique_ptr<Data> data = RcppUtilities::convert_data(train_matrix, sparse_train_matrix);

  data->set_outcome_index(outcome_index);
  data->set_treatment_index(treatment_index);
  data->set_instrument_index(treatment_index);
  if (use_sample_weights)
  {
    data->set_weight_index(sample_weight_index);
  }

  data->set_target_avg_weights(target_avg_weights);

  data->set_target_weight_penalty(target_weight_penalty);
  ...

}
```

During the tree splitting process, sample of input data and a list of variables from the input matrix X are  selected and passed to the function ``find_best_split`` . The function pre-computes and stores  the sum of treatment effect (``sum_node``) from all sample data. Note, the ``sample`` variable is a list of row index from the original data, and the function use it to  subset rows of  the input data . Then, the ``find_best_split_value`` finds the best splitting value from a variable of X by looping all variables ($x$) and all possible $x_{i}$ 

$$
sum\_node = \sum_{sample} sample\_weight_{i} * response\_by\_sample_{i} 
$$

```c++
// grf/src/splitting/InstrumentalSplittingRule.cpp.
bool InstrumentalSplittingRule::find_best_split(const Data &data,
                                                  size_t node,
                                                  const std::vector<size_t> &possible_split_vars,
                                                  const Eigen::ArrayXXd &responses_by_sample,
                                                  const std::vector<std::vector<size_t>> &samples,
                                                  std::vector<size_t> &split_vars,
                                                  std::vector<double> &split_values,
                                                  std::vector<bool> &send_missing_left)
  {
    // Precompute relevant quantities for this node.
    double weight_sum_node = 0.0;
    double sum_node = 0.0;
    double sum_node_z = 0.0;
    double sum_node_z_squared = 0.0;

    // calculate overall sum by add values of each sample
    for (auto &sample : samples[node])
    {
    //sample weights are given to each sample in estimation.
      double sample_weight = data.get_weight(sample); 
      weight_sum_node += sample_weight;

      //responses_by_sample is explained below
      sum_node += sample_weight * responses_by_sample(sample);

    // z is the treatment assignment (must be a binary or real numeric vector with no NAs).
      double z = data.get_instrument(sample);
      sum_node_z += sample_weight * z;
      sum_node_z_squared += sample_weight * z * z;
    }

    double size_node = sum_node_z_squared - sum_node_z * sum_node_z / weight_sum_node;
    double min_child_size = size_node * alpha;

    double mean_z_node = sum_node_z / weight_sum_node;
    size_t num_node_small_z = 0;
    for (auto &sample : samples[node])
    {
      double z = data.get_instrument(sample);
      if (z < mean_z_node)
      {
        num_node_small_z++;
      }
    }
    ...
  }
```

Here the ``responses_by_sample`` is a  pre-computed vector of all samples at the relabelling stage. 

$$
responses\_by\_sample_{i}   = (instrument_{i} - \bar{instrument})  * residual_{i} 
$$

$$
residual_{i}   = (response_{i} - \bar{response}) - local\_average\_treatment\_effect * (treatment_{i} - \bar{treatment}) 
$$

$$
local\_average\_treatment\_effect = \frac{\sum weight_{i} * (instrument_{i} - \bar{instrument}) * (outcome_{i} - \bar{outcome})}{\sum weight_{i} * (instrument_{i}  - \bar{instrument}) * (treatment_{i} - \bar{treatment})}
$$

```c++
//src/src/relabeling/InstrumentalRelabelingStrategy.cpp
bool InstrumentalRelabelingStrategy::relabel(
    const std::vector<size_t>& samples,
    const Data& data,
    Eigen::ArrayXXd& responses_by_sample) const {
    
    
    ...
    
      // Create the new outcomes per row from sample data
  for (size_t sample : samples) {
      // outcome is the Y input, treatment = instrument (binary variable)
    double response = data.get_outcome(sample);
    double treatment = data.get_treatment(sample);
    double instrument = data.get_instrument(sample);
    //reduced_form_weight is 0, so regularized_instrument = instrument
    double regularized_instrument = (1 - reduced_form_weight) * instrument + reduced_form_weight * treatment;
    // average_outcome is (sum of weight * outcome) / (sum of weight)
    //local_average_treatment_effect = sum ( weight * (instrument - average instrument) * (outcome  - average outcome) ) / sum(weight * (instrument -average instrument) * (treatment - average treatment))
    double residual = (response - average_outcome) - local_average_treatment_effect * (treatment - average_treatment);

    //average_regularized_instrument is sample average of instrument 
    responses_by_sample(sample) = (regularized_instrument - average_regularized_instrument) * residual;
  }
}
```

The process of finding the best splitting point is defined below
```c++
// grf/src/splitting/InstrumentalSplittingRule.cpp.
 void InstrumentalSplittingRule::find_best_split_value(const Data &data,
                                                        size_t node, size_t var,
                                                        size_t num_samples,
                                                        double weight_sum_node,
                                                        double sum_node,
                                                        double mean_node_z,
                                                        size_t num_node_small_z,
                                                        double sum_node_z,
                                                        double sum_node_z_squared,
                                                        double min_child_size,
                                                        double &best_value,
                                                        size_t &best_var,
                                                        double &best_decrease,
                                                        bool &best_send_missing_left,
                                                        const Eigen::ArrayXXd &responses_by_sample,
                                                        const std::vector<std::vector<size_t>> &samples)
  
    {...}
```

It first pre-computes and stores total treatment effect by looping over all value in the variable. If the current value is different from the next value, then it adds one more splitting index.

The goal is to find the best splitting point ($i$) of  $x_{j}$ that maximize $decrease_{i}$:

$$
decrease_{i} = \frac{ sum\_left_{i}^{2}}{weight\_sum\_left_{i}}  + \frac{ sum\_right_{i}^{2}}{weight\_sum\_right_{i}} 
$$

sum_left is the sum of all weighted response ( $j$ ) on the samples in left child node. 

$$
\begin{aligned}
    sum\_left_{i} &= \sum_{j} sample\_weight_{j} * response\_by\_sample_{j} \\
    weight\_sum\_left_{i} &= \sum_{j} sample\_weight_{j}
\end{aligned}
$$

Then the right values are calculated by total values - left values:

$$
\begin{aligned}
    sum\_left_{i} &= sum\_node - sum\_left \\
    weight\_sum\_right_{i} &= weight\_sum\_node - weight\_sum\_left 
\end{aligned}
$$

where sum_node is the total value from all sample data ($j$) before finding best splitting point.

$$
sum\_node = \sum_{j} sample\_weight_{j} * response\_by\_sample_{j}
$$

```c++
    // var is an index column, 0 is the first variable of input X
    // sorted_samples is an ordered vector for selected variable x 
    size_t split_index = 0;
    for (size_t i = 0; i < num_samples - 1; i++)
    {
        size_t sample = sorted_samples[i];
        size_t next_sample = sorted_samples[i + 1];
        double sample_value = data.get(sample, var);
        double z = data.get_instrument(sample);
        double sample_weight = data.get_weight(sample);

        weight_sums[split_index] += sample_weight; //vector
        sums[split_index] += sample_weight * responses_by_sample(sample);//vector
        ++counter[split_index];//vector

        sums_z[split_index] += sample_weight * z;
        sums_z_squared[split_index] += sample_weight * z * z;
        if (z < mean_node_z)
        {
            ++num_small_z[split_index]; //vector
        }
        // if the next sample value is different, then add splitting index 
        if (sample_value != next_sample_value && !std::isnan(next_sample_value))
        {
            ++split_index;
        }
  }

```

For each splitting index, compute the decrease of impurity on both left and right child node. We only need to calculate the left treatment effects, then we can get the right values by subtracting the left from the total value.

```c++

     for (size_t i = 0; i < num_splits; ++i)
      {
        n_left += counter[i]; // number of values in the left
        num_left_small_z += num_small_z[i];
        weight_sum_left += weight_sums[i];
        sum_left += sums[i];
        sum_left_z += sums_z[i];
        sum_left_z_squared += sums_z_squared[i];

        // Skip this split if the left child does not contain enough
        // z values below and above the parent's mean.
        size_t num_left_large_z = n_left - num_left_small_z;
        if (num_left_small_z < min_node_size || num_left_large_z < min_node_size)
        {
          continue;
        }

        // Stop if the right child does not contain enough z values below
        // and above the parent's mean.
        size_t n_right = num_samples - n_left;
        size_t num_right_small_z = num_node_small_z - num_left_small_z;
        size_t num_right_large_z = n_right - num_right_small_z;
        if (num_right_small_z < min_node_size || num_right_large_z < min_node_size)
        {
          break;
        }

        // Calculate relevant quantities for the left child.
        double size_left = sum_left_z_squared - sum_left_z * sum_left_z / weight_sum_left;
        // Skip this split if the left child's variance is too small.
        if (size_left < min_child_size || (imbalance_penalty > 0.0 && size_left == 0))
        {
          continue;
        }

        // Calculate relevant quantities for the right child.
        double weight_sum_right = weight_sum_node - weight_sum_left;
        double sum_right = sum_node - sum_left;
        double sum_right_z_squared = sum_node_z_squared - sum_left_z_squared;
        double sum_right_z = sum_node_z - sum_left_z;
        double size_right = sum_right_z_squared - sum_right_z * sum_right_z / weight_sum_right;

        // Skip this split if the right child's variance is too small.
        if (size_right < min_child_size || (imbalance_penalty > 0.0 && size_right == 0))
        {
          continue;
        }

        // Calculate the decrease in impurity.
        double decrease_left = sum_left * sum_left / weight_sum_left;
        double decrease_right = sum_right * sum_right / weight_sum_right;
        double decrease = decrease_left + decrease_right;
        // Penalize splits that are too close to the edges of the data.
        decrease -= imbalance_penalty * (1.0 / size_left + 1.0 / size_right);
        if (target_weight_penalty > 0)
        { // my contribution (explained in the next section)
          Eigen::VectorXd left_target_avg_weight = target_avg_weights_sorted.topRows(n_left).colwise().mean();
          Eigen::VectorXd right_target_avg_weight = target_avg_weights_sorted.bottomRows(n_right).colwise().mean();

          double penalty_target_weight = left_target_avg_weight.lpNorm<2>() * decrease_left + right_target_avg_weight.lpNorm<2>() * decrease_right; /// weight_sum_left   / weight_sum_right
          // std::cout << "var" << var << "decrease:" << decrease << "penalty:" << penalty_target_weight << "\n";
          decrease -= penalty_target_weight * target_weight_penalty;
        }

  
        // store the current best value and splitting point , default is 0
        if (decrease > best_decrease)
        {

          best_value = possible_split_values[i];
          best_var = var;
          best_decrease = decrease;
          best_send_missing_left = send_left;
        }
}
```


## Target Weight Penalty 
