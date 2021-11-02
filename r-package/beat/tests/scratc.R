library(Rcpp)
sourceCpp('test.cpp')


a = matrix(seq(1:120), c(20, 4))
b = list(a,a,a)

get_matrix(b)


get_length(b)
get_vector(b, 0, 1)

conv(a[1,] )

get_ncols(a, 2)
 head(a, 2)
 cppFunction("bool conv(NumericVector X) { \
      Eigen::Map<Eigen::VectorXd>  XS(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(X)); return true; } ",
             depends="RcppEigen")
 conv(1:4)
