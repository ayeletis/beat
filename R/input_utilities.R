#' @import data.table


validate_X <- function(X, allow.na = FALSE) {
    valid.classes <- c("matrix", "data.frame", "dgCMatrix")

    if (!inherits(X, valid.classes)) {
        stop(paste(
      "Currently the only supported data input types are:",
      "`matrix`, `data.frame`, `dgCMatrix`"
    ))
    }

    if (!inherits(X, "dgCMatrix") && !is.numeric(as.matrix(X))) {
        stop(paste(
      "The feature matrix X must be numeric. GRF does not",
      "currently support non-numeric features. If factor variables",
      "are required, we recommend one of the following: Either",
      "represent the factor with a 1-vs-all expansion,",
      "(e.g., using model.matrix(~. , data=X)), or then encode the factor",
      "as a numeric via any natural ordering (e.g., if the factor is a month).",
      "For more on GRF and categorical variables see the online vignette:",
      "https://grf-labs.github.io/grf/articles/categorical_inputs.html"
    ))
    }

    has.missing.values <- any(is.na(X))

    if (!allow.na && has.missing.values) {
        stop("The feature matrix X contains at least one NA.")
    }

    has.missing.values
}

validate_observations <- function(V, X, allow.matrix = FALSE) {
    if (!allow.matrix) {
        if (is.matrix(V) && ncol(V) == 1) {
            V <- as.vector(V)
        } else if (!is.vector(V)) {
            stop(paste("Observations (W, Y, Z or D) must be vectors."))
        }
    }

    if (!is.numeric(V) && !is.logical(V)) {
        stop(paste(
      "Observations (W, Y, Z or D) must be numeric. GRF does not ",
      "currently support non-numeric observations."
    ))
    }

    if (any(is.na(V))) {
        stop("The vector of observations (W, Y, Z or D) contains at least one NA.")
    }

    if (NROW(V) != nrow(X)) {
        stop("length of observation (W, Y, Z or D) does not equal nrow(X).")
    }
    V
}

validate_num_threads <- function(num.threads) {
    if (is.null(num.threads)) {
        num.threads <- 0
    } else if (!is.numeric(num.threads) | num.threads < 0) {
        stop("Error: Invalid value for num.threads")
    }
    num.threads
}

validate_clusters <- function(clusters, X) {
    if (is.null(clusters) || length(clusters) == 0) {
        return(vector(mode = "numeric", length = 0))
    }
    if (mode(clusters) != "numeric") {
        stop("Clusters must be able to be coerced to a numeric vector.")
    }
    clusters <- as.numeric(clusters)
    if (!all(clusters == floor(clusters))) {
        stop("Clusters vector cannot contain floating point values.")
    } else if (length(clusters) != nrow(X)) {
        stop("Clusters vector has incorrect length.")
    } else {
        # convert to integers between 0 and n clusters
        clusters <- as.numeric(as.factor(clusters)) - 1
    }
    clusters
}

validate_equalize_cluster_weights <- function(equalize.cluster.weights, clusters, sample.weights) {
    if (is.null(clusters) || length(clusters) == 0) {
        return(0)
    }
    cluster_size_counts <- table(clusters)
    if (equalize.cluster.weights == TRUE) {
        samples.per.cluster <- min(cluster_size_counts)
        if (!is.null(sample.weights)) {
            stop("If equalize.cluster.weights is TRUE, sample.weights must be NULL.")
        }
    } else if (equalize.cluster.weights == FALSE) {
        samples.per.cluster <- max(cluster_size_counts)
    } else {
        stop("equalize.cluster.weights must be either TRUE or FALSE.")
    }

    samples.per.cluster
}

validate_boost_error_reduction <- function(boost.error.reduction) {
    if (boost.error.reduction < 0 || boost.error.reduction > 1) {
        stop("boost.error.reduction must be between 0 and 1")
    }
    boost.error.reduction
}

validate_ll_vars <- function(linear.correction.variables, num.cols) {
    if (is.null(linear.correction.variables)) {
        linear.correction.variables <- 1:num.cols
    }
    if (min(linear.correction.variables) < 1) {
        stop("Linear correction variables must take positive integer values.")
    } else if (max(linear.correction.variables) > num.cols) {
        stop("Invalid range of correction variables.")
    } else if (!is.vector(linear.correction.variables) ||
    !all(linear.correction.variables == floor(linear.correction.variables))) {
        stop("Linear correction variables must be a vector of integers.")
    }
    linear.correction.variables
}

validate_ll_lambda <- function(lambda) {
    if (lambda < 0) {
        stop("Lambda cannot be negative.")
    } else if (!is.numeric(lambda) || length(lambda) > 1) {
        stop("Lambda must be a scalar.")
    }
    lambda
}

validate_ll_path <- function(lambda.path) {
    if (is.null(lambda.path)) {
        lambda.path <- c(0, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 1, 10)
    } else if (min(lambda.path) < 0) {
        stop("Lambda values cannot be negative.")
    } else if (!is.numeric(lambda.path)) {
        stop("Lambda values must be numeric.")
    }
    lambda.path
}

validate_ll_cutoff <- function(ll.split.cutoff, num.rows) {
    if (is.null(ll.split.cutoff)) {
        ll.split.cutoff <- floor(sqrt(num.rows))
    } else if (!is.numeric(ll.split.cutoff) || length(ll.split.cutoff) > 1) {
        stop("LL split cutoff must be NULL or a scalar")
    } else if (ll.split.cutoff < 0 || ll.split.cutoff > num.rows) {
        stop("Invalid range for LL split cutoff")
    }
    ll.split.cutoff
}

validate_newdata <- function(newdata, X, allow.na = FALSE) {
    validate_X(newdata, allow.na = allow.na)
    if (ncol(newdata) != ncol(X)) {
        stop("newdata must have the same number of columns as the training matrix.")
    }
}

validate_sample_weights <- function(sample.weights, X) {
    if (!is.null(sample.weights)) {
        if (length(sample.weights) != nrow(X)) {
            stop("sample.weights has incorrect length")
        }
        if (any(sample.weights < 0)) {
            stop("sample.weights must be nonnegative")
        }
    }
}

#' @importFrom Matrix Matrix cBind
#' @importFrom methods new
# Indices are offset by 1 for C++.
create_train_matrices <- function(X,
                                  outcome = NULL,
                                  treatment = NULL,
                                  instrument = NULL,
                                  survival.numerator = NULL,
                                  survival.denominator = NULL,
                                  censor = NULL,
                                  sample.weights = FALSE) {
    default.data <- matrix(nrow = 0, ncol = 0)
    sparse.data <- new("dgCMatrix", Dim = c(0L, 0L))
    out <- list()
    offset <- ncol(X) - 1
    if (!is.null(outcome)) {
        out[["outcome.index"]] <- (offset + 1):(offset + NCOL(outcome))
        offset <- offset + NCOL(outcome)
    }
    if (!is.null(treatment)) {
        out[["treatment.index"]] <- offset + 1
        offset <- offset + 1
    }
    if (!is.null(instrument)) {
        out[["instrument.index"]] <- offset + 1
        offset <- offset + 1
    }
    if (!is.null(survival.numerator)) {
        out[["causal.survival.numerator.index"]] <- offset + 1
        offset <- offset + 1
    }
    if (!is.null(survival.denominator)) {
        out[["causal.survival.denominator.index"]] <- offset + 1
        offset <- offset + 1
    }
    if (!is.null(censor)) {
        out[["censor.index"]] <- offset + 1
        offset <- offset + 1
    }
    # Forest bindings without sample weights: sample.weights = FALSE
    # Forest bindings with sample weights:
    # -sample.weights = NULL if no weights passed
    # -sample.weights = numeric vector if passed
    if (is.logical(sample.weights)) {
        sample.weights <- NULL
    } else {
        out[["sample.weight.index"]] <- offset + 1
        if (is.null(sample.weights)) {
            out[["use.sample.weights"]] <- FALSE
        } else {
            out[["use.sample.weights"]] <- TRUE
        }
    }

    if (inherits(X, "dgCMatrix") && ncol(X) > 1) {
        sparse.data <- cbind(X, outcome, treatment, instrument, survival.numerator, survival.denominator, censor, sample.weights)
    } else {
        X <- as.matrix(X)
        default.data <- as.matrix(cbind(X, outcome, treatment, instrument, survival.numerator, survival.denominator, censor, sample.weights))
    }
    out[["train.matrix"]] <- default.data
    out[["sparse.train.matrix"]] <- sparse.data

    out
}

#' @importFrom Matrix Matrix cBind
#' @importFrom methods new
create_test_matrices <- function(X) {
    default.data <- matrix(nrow = 0, ncol = 0)
    sparse.data <- new("dgCMatrix", Dim = c(0L, 0L))
    out <- list()
    if (inherits(X, "dgCMatrix") && ncol(X) > 1) {
        sparse.data <- X
    } else {
        default.data <- as.matrix(X)
    }
    out[["test.matrix"]] <- default.data
    out[["sparse.test.matrix"]] <- sparse.data

    out
}

observation_weights <- function(forest) {
    # Case 1: No sample.weights
    if (is.null(forest$sample.weights)) {
        if (length(forest$clusters) == 0 || !forest$equalize.cluster.weights) {
            raw.weights <- rep(1, length(forest$Y.orig))
        } else {
            # If clustering with no sample.weights provided and equalize.cluster.weights = TRUE, then
            # give each observation weight 1/cluster size, so that the total weight of each cluster is the same.
            clust.factor <- factor(forest$clusters)
            inverse.counts <- 1 / as.numeric(Matrix::colSums(Matrix::sparse.model.matrix(~clust.factor + 0)))
            raw.weights <- inverse.counts[as.numeric(clust.factor)]
        }
    }

    # Case 2: sample.weights provided
    if (!is.null(forest$sample.weights)) {
        if (length(forest$clusters) == 0 || !forest$equalize.cluster.weights) {
            raw.weights <- forest$sample.weights
        } else {
            stop("Specifying non-null sample.weights is not allowed when equalize.cluster.weights = TRUE")
        }
    }

    return(raw.weights / sum(raw.weights))
}

# Call the grf Rcpp bindings (argument_names) with R argument.names
#
# All the bindings argument names (C++) have underscores: sample_weights, train_matrix, etc.
# On the R side each variable name is written as sample.weights, train.matrix, etc.
# This function simply replaces the underscores in the passed argument names with dots.
do.call.rcpp = function(what, args, quote = FALSE, envir = parent.frame()) {
    names(args) = gsub("\\.", "_", names(args))
    do.call(what, args, quote, envir)
}

validate_subset <- function(forest, subset) {
    if (is.null(subset)) {
        subset <- 1:length(forest$Y.orig)
    }
    if (class(subset) == "logical" && length(subset) == length(forest$Y.orig)) {
        subset <- which(subset)
    }
    if (!all(subset %in% 1:length(forest$Y.orig))) {
        stop(paste(
      "If specified, subset must be a vector contained in 1:n,",
      "or a boolean vector of length n."
    ))
    }
    subset
}


standardize = function(x) {
    if (length(unique(x)) == 1) {
        return(scale(x, center = TRUE, scale = FALSE))
    } else {
        return(scale(x, center = TRUE, scale = TRUE))
    }
}

#'
construct_target_weight_mean = function(x, z, num_breaks = 256) {
    calculate_avg = function(dat, x_col, z_col, num_breaks) {
        df = dat[, .SD, .SDcols = c(x_col, z_col)]
        df[, cut_bins := cut(get(x_col), breaks = num_breaks)]

        df[, mean_value := mean(get(z_col)), by = cut_bins]
        out = df[, mean_value]

        return(out)
    }

    stopifnot(is.matrix(z))
    stopifnot(dim(x)[1] == dim(z)[1])
    dat = as.data.table(x)
    x_cols = paste0('x_', 1:dim(x)[2])
    names(dat) = x_cols
    dat_z = as.data.table(z)
    z_cols = paste0("z_", 1:dim(z)[2])
    names(dat_z) = z_cols
    dat = cbind(dat, dat_z)
    # output matrix: 3D [var, n, z]
    out = vector('list', length = length(x_cols))
    for (i in 1:length(x_cols)) {
        #  [# observation, # z columns], sapply combine by columns
        x_z = sapply(z_cols, calculate_avg, x_col = x_cols[i], num_breaks = num_breaks, dat = dat)
        out[[i]] = t(x_z) # arma::mat is column base
    }

    target.avg.weights = array(unlist(out), dim = c(dim(out[[1]]), length(out)))

    return(target.avg.weights)
}

#'
available_metrics = c("split_l2_norm_rate", # left, right: l2 norm(colmean target weight)* penalty rate * node decrease
                        "euclidean_distance_rate", # (left+right decrease) *  Euclidean distance (column mean target weight left, right ) * penalty rate
                        "cosine_similarity_rate", # (left+right decrease) *  (1-cos_sim(column mean target weight left, right )) * penalty rate

                        "split_l2_norm", #  sum(left,right l2 norm(colmean target weight))* penalty rate
                        "euclidean_distance", #  Euclidean distance (column mean target weight left, right ) * penalty rate
                        "cosine_similarity" #  (1-cos_sim(column mean target weight left, right )) * penalty rate
                        )

#'
get_available_distance_metrics = function() {
    return(available_metrics)
}