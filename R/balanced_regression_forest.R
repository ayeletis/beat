#' Balanced Regression forest
#'
#' Trains a regression forest that can be used to estimate
#' the conditional mean function mu(x) = E[Y | X = x]
#'
#' @param X The covariates used in the regression.
#' @param Y The outcome.
#' @param num.trees Number of trees grown in the forest. Note: Getting accurate
#'                  confidence intervals generally requires more trees than
#'                  getting accurate predictions. Default is 2000.
#' @param sample.weights Weights given to an observation in estimation.
#'                       If NULL, each observation is given the same weight. Default is NULL.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'  Default is NULL (ignored).
#' @param equalize.cluster.weights If FALSE, each unit is given the same weight (so that bigger
#'  clusters get more weight). If TRUE, each cluster is given equal weight in the forest. In this case,
#'  during training, each tree uses the same number of observations from each drawn cluster: If the
#'  smallest cluster has K units, then when we sample a cluster during training, we only give a random
#'  K elements of the cluster to the tree-growing procedure. When estimating average treatment effects,
#'  each observation is given weight 1/cluster size, so that the total weight of each cluster is the
#'  same. Note that, if this argument is FALSE, sample weights may also be directly adjusted via the
#'  sample.weights argument. If this argument is TRUE, sample.weights must be set to NULL. Default is
#'  FALSE.
#' @param sample.fraction Fraction of the data used to build each tree.
#'                        Note: If honesty = TRUE, these subsamples will
#'                        further be cut by a factor of honesty.fraction. Default is 0.5.
#' @param mtry Number of variables tried for each split. Default is
#'             \eqn{\sqrt p + 20} where p is the number of variables.
#' @param min.node.size A target for the minimum number of observations in each tree leaf. Note that nodes
#'                      with size smaller than min.node.size can occur, as in the original randomForest package.
#'                      Default is 5.
#' @param honesty Whether to use honest splitting (i.e., sub-sample splitting). Default is TRUE.
#'  For a detailed description of honesty, honesty.fraction, honesty.prune.leaves, and recommendations for
#'  parameter tuning, see the grf
#'  \href{https://grf-labs.github.io/grf/REFERENCE.html#honesty-honesty-fraction-honesty-prune-leaves}{algorithm reference}.
#' @param honesty.fraction The fraction of data that will be used for determining splits if honesty = TRUE. Corresponds
#'                         to set J1 in the notation of the paper. Default is 0.5 (i.e. half of the data is used for
#'                         determining splits).
#' @param honesty.prune.leaves If TRUE, prunes the estimation sample tree such that no leaves
#'  are empty. If FALSE, keep the same tree as determined in the splits sample (if an empty leave is encountered, that
#'  tree is skipped and does not contribute to the estimate). Setting this to FALSE may improve performance on
#'  small/marginally powered data, but requires more trees (note: tuning does not adjust the number of trees).
#'  Only applies if honesty is enabled. Default is TRUE.
#' @param alpha A tuning parameter that controls the maximum imbalance of a split. Default is 0.05.
#' @param imbalance.penalty A tuning parameter that controls how harshly imbalanced splits are penalized. Default is 0.
#' @param ci.group.size The forest will grow ci.group.size trees on each subsample.
#'                      In order to provide confidence intervals, ci.group.size must
#'                      be at least 2. Default is 2.
#' @param tune.parameters A vector of parameter names to tune.
#'  If "all": all tunable parameters are tuned by cross-validation. The following parameters are
#'  tunable: ("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
#'   "honesty.prune.leaves", "alpha", "imbalance.penalty"). If honesty is FALSE the honesty.* parameters are not tuned.
#'  Default is "none" (no parameters are tuned).
#' @param tune.num.trees The number of trees in each 'mini forest' used to fit the tuning model. Default is 50.
#' @param tune.num.reps The number of forests used to fit the tuning model. Default is 100.
#' @param tune.num.draws The number of random parameter values considered when using the model
#'                          to select the optimal parameters. Default is 1000.
#' @param compute.oob.predictions Whether OOB predictions on training set should be precomputed. Default is TRUE.
#' @param num.threads Number of threads used in training. By default, the number of threads is set
#'                    to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A trained regression forest object. If tune.parameters is enabled,
#'  then tuning information will be included through the `tuning.output` attribute.
#'
#' @examples
#' \donttest{
#' # Train a standard regression forest.
#' n <- 500
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' r.forest <- regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' r.pred <- predict(r.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred <- predict(r.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' r.forest <- regression_forest(X, Y, num.trees = 100)
#' r.pred <- predict(r.forest, X.test, estimate.variance = TRUE)
#' }
#'
#'
#' @import Rcpp
#'
#' @export
#' @importFrom utils modifyList
balanced_regression_forest <- function(X, Y,
                              num.trees = 2000,
                              sample.weights = NULL,
                              target.weights = NULL,
                              target.weight.penalty = 0,
                              target.weight.bins.breaks = 256,
                              target.weight.standardize = TRUE,
                              target.avg.weights = NULL,
                              clusters = NULL,
                              equalize.cluster.weights = FALSE,
                              sample.fraction = 0.5,
                              mtry = min(ceiling(sqrt(ncol(X)) + 20), ncol(X)),
                              min.node.size = 5,
                              honesty = TRUE,
                              honesty.fraction = 0.5,
                              honesty.prune.leaves = TRUE,
                              alpha = 0.05,
                              imbalance.penalty = 0,
                              ci.group.size = 2,
                              tune.parameters = "none",
                              tune.num.trees = 50,
                              tune.num.reps = 100,
                              tune.num.draws = 1000,
                              compute.oob.predictions = TRUE,
                              num.threads = NULL,
                              seed = runif(1, 0, .Machine$integer.max)) {
  # target.avg.weights: array, [num target weight, observation, num x]
  # if null, create it internally
    has.missing.values <- validate_X(X, allow.na = TRUE)
    validate_sample_weights(sample.weights, X)
    Y <- validate_observations(Y, X)
    stopifnot(target.weight.penalty >= 0)
    stopifnot(nrow(target.weights) == nrow(X))

    # target weight
    if (!is.null(target.weights)) {
        stopifnot(dim(target.weights)[1] == dim(X)[1])
        if (isTRUE(target.weight.standardize)) {
            target.weights = apply(target.weights, 2, standardize) # per column
        }

    } else {
        #target.weights = as.matrix(replicate(NROW(X), 0))
        stop("If target.weight is missing, use regression_forest")

    }
    # verify penalty metric
    available_metrics = c("split_l2_norm_rate", # left, right: l2 norm(colmean target weight)* penalty rate * node decrease
                        "euclidean_distance_rate", # (left+right decrease) *  Euclidean distance (column mean target weight left, right ) * penalty rate
                        "cosine_similarity_rate", # (left+right decrease) *  (1-cos_sim(column mean target weight left, right )) * penalty rate

                        "split_l2_norm", #  sum(left,right l2 norm(colmean target weight))* penalty rate
                        "euclidean_distance", #  Euclidean distance (column mean target weight left, right ) * penalty rate
                        "cosine_similarity" #  (1-cos_sim(column mean target weight left, right )) * penalty rate
                        )

    if(is.null(target.avg.weights)){
      target.avg.weights = construct_target_weight_mean(x = X, z = target.weights,
                                                        num_breaks = target.weight.bins.breaks)
    }
    stopifnot(is.array(target.avg.weights))
    stopifnot(dim(target.avg.weights) == c(dim(target.weights)[2], dim(X)[1], dim(X)[2]))

    #output list : [dim(X)[2]] [[num target weight feature, num rows obs]]
    # then  convert to 3d array: [dim(target weight), length(list)]

    ##


    clusters <- validate_clusters(clusters, X)
    samples.per.cluster <- validate_equalize_cluster_weights(equalize.cluster.weights, clusters, sample.weights)
    num.threads <- validate_num_threads(num.threads)

    all.tunable.params <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
                          "honesty.prune.leaves", "alpha", "imbalance.penalty", "target.weight.penalty")

    if (max(abs(target.weights)) == 0) {
        all.tunable.params = all.tunable.params[all.tunable.params != 'target.weight.penalty']
    }

    data <- create_train_matrices(X, outcome = Y, sample.weights = sample.weights)
    args <- list(num.trees = num.trees,
              target.avg.weights = target.avg.weights,
              target.weight.penalty = target.weight.penalty,
              target.weight.penalty.metric = "split_l2_norm_rate",
               clusters = clusters,
               samples.per.cluster = samples.per.cluster,
               sample.fraction = sample.fraction,
               mtry = mtry,
               min.node.size = min.node.size,
               honesty = honesty,
               honesty.fraction = honesty.fraction,
               honesty.prune.leaves = honesty.prune.leaves,
               alpha = alpha,
               imbalance.penalty = imbalance.penalty,
               ci.group.size = ci.group.size,
               compute.oob.predictions = compute.oob.predictions,
               num.threads = num.threads,
               seed = seed)

    tuning.output <- NULL
    if (!identical(tune.parameters, "none")) {
        tuning.output <- tune_balanced_regression_forest(X, Y,
                                            sample.weights = sample.weights,
                                            target.avg.weights = target.avg.weights,
                                            target.weight.penalty = target.weight.penalty,
                                            target.weight.penalty.metric = "split_l2_norm_rate",
                                            clusters = clusters,
                                            equalize.cluster.weights = equalize.cluster.weights,
                                            sample.fraction = sample.fraction,
                                            mtry = mtry,
                                            min.node.size = min.node.size,
                                            honesty = honesty,
                                            honesty.fraction = honesty.fraction,
                                            honesty.prune.leaves = honesty.prune.leaves,
                                            alpha = alpha,
                                            imbalance.penalty = imbalance.penalty,
                                            ci.group.size = ci.group.size,
                                            tune.parameters = tune.parameters,
                                            tune.num.trees = tune.num.trees,
                                            tune.num.reps = tune.num.reps,
                                            tune.num.draws = tune.num.draws,
                                            num.threads = num.threads,
                                            seed = seed)
        args <- modifyList(args, as.list(tuning.output[["params"]]))
    }

    forest <- do.call.rcpp(balanced_regression_train, c(data, args))
    class(forest) <- c("balanced_regression_forest", "regression_forest", "grf")
    forest[["ci.group.size"]] <- ci.group.size
    forest[["X.orig"]] <- X
    forest[["Y.orig"]] <- Y
    forest[["sample.weights"]] <- sample.weights
    forest[["clusters"]] <- clusters
    forest[["equalize.cluster.weights"]] <- equalize.cluster.weights
    forest[["tunable.params"]] <- args[all.tunable.params]
    forest[["tuning.output"]] <- tuning.output
    forest[["has.missing.values"]] <- has.missing.values
    forest[["target.weight.penalty.metric"]] <- "split_l2_norm_rate"

    forest
}

#' Predict with a regression forest
#'
#' Gets estimates of E[Y|X=x] using a trained regression forest.
#'
#' @param object The trained forest.
#' @param newdata Points at which predictions should be made. If NULL, makes out-of-bag
#'                predictions on the training set instead (i.e., provides predictions at
#'                Xi using only trees that did not use the i-th training example). Note
#'                that this matrix should have the number of columns as the training
#'                matrix, and that the columns must appear in the same order.
#' @param linear.correction.variables Optional subset of indexes for variables to be used in local
#'                   linear prediction. If NULL, standard GRF prediction is used. Otherwise,
#'                   we run a locally weighted linear regression on the included variables.
#'                   Please note that this is a beta feature still in development, and may slow down
#'                   prediction considerably. Defaults to NULL.
#' @param ll.lambda Ridge penalty for local linear predictions
#' @param ll.weight.penalty Option to standardize ridge penalty by covariance (TRUE),
#'                            or penalize all covariates equally (FALSE). Defaults to FALSE.
#' @param num.threads Number of threads used in training. If set to NULL, the software
#'                    automatically selects an appropriate amount.
#' @param estimate.variance Whether variance estimates for hat{tau}(x) are desired
#'                          (for confidence intervals).
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of predictions, along with estimates of the error and
#'         (optionally) its variance estimates. Column 'predictions' contains
#'         estimates of E[Y|X=x]. The square-root of column 'variance.estimates' is the standard error
#          of these predictions. Column 'debiased.error' contains out-of-bag estimates of
#'         the test mean-squared error. Column 'excess.error' contains
#'         jackknife estimates of the Monte-carlo error. The sum of 'debiased.error'
#'         and 'excess.error' is the raw error attained by the current forest, and
#'         'debiased.error' alone is an estimate of the error attained by a forest with
#'         an infinite number of trees. We recommend that users grow
#'         enough forests to make the 'excess.error' negligible.
#'
#' @examples
#' \donttest{
#' # Train a standard regression forest.
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1] * rnorm(n)
#' r.forest <- regression_forest(X, Y)
#'
#' # Predict using the forest.
#' X.test <- matrix(0, 101, p)
#' X.test[, 1] <- seq(-2, 2, length.out = 101)
#' r.pred <- predict(r.forest, X.test)
#'
#' # Predict on out-of-bag training samples.
#' r.pred <- predict(r.forest)
#'
#' # Predict with confidence intervals; growing more trees is now recommended.
#' r.forest <- regression_forest(X, Y, num.trees = 100)
#' r.pred <- predict(r.forest, X.test, estimate.variance = TRUE)
#' }
#'
#' @method predict balanced_regression_forest
#' @export
predict.balanced_regression_forest <- function(object, newdata = NULL,
                                      linear.correction.variables = NULL,
                                      ll.lambda = NULL,
                                      ll.weight.penalty = FALSE,
                                      num.threads = NULL,
                                      estimate.variance = FALSE,
                                      ...) {
    local.linear <- !is.null(linear.correction.variables)
    allow.na <- !local.linear

    # If possible, use pre-computed predictions.
    if (is.null(newdata) && !estimate.variance && !local.linear && !is.null(object$predictions)) {
        return(data.frame(
      predictions = object$predictions,
      debiased.error = object$debiased.error,
      excess.error = object$excess.error))
    }

    num.threads <- validate_num_threads(num.threads)
    forest.short <- object[-which(names(object) == "X.orig")]
    X <- object[["X.orig"]]
    train.data <- create_train_matrices(X, outcome = object[["Y.orig"]])

    if (local.linear) {
        linear.correction.variables <- validate_ll_vars(linear.correction.variables, ncol(X))

        if (is.null(ll.lambda)) {
            ll.regularization.path <- tune_ll_regression_forest(
        object, linear.correction.variables,
        ll.weight.penalty, num.threads)
            ll.lambda <- ll.regularization.path$lambda.min
        } else {
            ll.lambda <- validate_ll_lambda(ll.lambda)
        }

        # subtract 1 to account for C++ indexing
        linear.correction.variables <- linear.correction.variables - 1
    }
    args <- list(forest.object = forest.short,
               num.threads = num.threads,
               estimate.variance = estimate.variance)
    ll.args <- list(ll.lambda = ll.lambda,
                  ll.weight.penalty = ll.weight.penalty,
                  linear.correction.variables = linear.correction.variables)

    if (!is.null(newdata)) {
        validate_newdata(newdata, X, allow.na = allow.na)
        test.data <- create_test_matrices(newdata)
        if (!local.linear) {
            ret <- do.call.rcpp(regression_predict, c(train.data, test.data, args))
        } else {
            ret <- do.call.rcpp(ll_regression_predict, c(train.data, test.data, args, ll.args))
        }
    } else {
        if (!local.linear) {
            ret <- do.call.rcpp(regression_predict_oob, c(train.data, args))
        } else {
            ret <- do.call.rcpp(ll_regression_predict_oob, c(train.data, args, ll.args))
        }
    }

    # Convert list to data frame.
    empty <- sapply(ret, function(elem) length(elem) == 0)
    do.call(cbind.data.frame, ret[!empty])
}
