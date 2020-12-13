#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
std::vector<Eigen::MatrixXd> get_matrix(Rcpp::List x){
  std::vector<Eigen::MatrixXd> M(x.size());
  for(int i; i<x.size(); i++){
    NumericMatrix a = x[i];
    Eigen::MatrixXd Xs = Rcpp::as<Eigen::MatrixXd>(a);
    M[i] = Xs;
  }
  return M;
}



// [[Rcpp::export]]
Eigen::VectorXf conv(NumericMatrix  X, size_t n ){
  NumericVector a = X.row(n);
  Eigen::VectorXf Xs = Rcpp::as<Eigen::VectorXf>(a);

  // Eigen::Map<Eigen::VectorXf> XS(Rcpp::as<Eigen::Map<Eigen::VectorXf> >(X));
  return Xs*2;
}


// [[Rcpp::export]]
Eigen::MatrixXf get_ncols(NumericMatrix x, size_t n){
  NumericMatrix M(n, x.ncol());
  for(size_t t=0; t <n; t++){
    M.row(t) = x.row(t);
  }
  Eigen::MatrixXf Y(Rcpp::as<Eigen::MatrixXf>(M));

  // for(size_t j=0; j<n; j++){
  //   NumericVector a = M.row(j)
  //   Eigen::VectorXf Xs = Rcpp::as<Eigen::VectorXf>(a);
  //   Y.row(j) =  Xs;
  // }

  return Y;
}

// [[Rcpp::export]]
Rcpp::NumericVector get_vector(Rcpp::List x, size_t n, size_t row){
  Rcpp::NumericMatrix M = x[n];
  return M.row(row);
}


// [[Rcpp::export]]
std::vector<Eigen::MatrixXf>  get_one_row(Rcpp::List x){
  // Rcpp::NumericMatrix M = x[0];
  // const Eigen::MatrixXf Xs = Rcpp::as<Eigen::MatrixXf>(M);

  std::vector<Eigen::MatrixXf> Y;
  std::vector<Eigen::MatrixXf> OUT(20);
  // OUT[0] = Xs;
  for(size_t i; i<x.size(); i++){
    Rcpp::NumericMatrix M = x[i];
    // Eigen::MatrixXf Xs = Rcpp::as<Eigen::MatrixXf>(M);
    OUT[i] = Rcpp::as<Eigen::MatrixXf>(M);//Xs;
  }
  Y = OUT;
  return Y;
}


// [[Rcpp::export]]
void check_values(Rcpp::NumericMatrix x){
  Rcpp::Rcout<< x.rows() << "\n";
  Rcpp::Rcout<< x.cols() << "\n";

  for(double i; i<x.rows(); i++){
    for(double j; j<x.cols(); j++){
      Rcpp::Rcout << "Value: " << x(i,j) << "\n";
      Rcout << "The value is \n";
    }
  }
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

