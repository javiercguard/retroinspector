#include <Rcpp.h>
using namespace Rcpp;

double binomSniffles (double x, double size, double chance) {
  return pow(chance, x) * pow((1 - chance), (size - x));
}

// [[Rcpp::export]]
NumericVector dbinomSnifflesCpp(double x, double size, NumericVector chance) {
  NumericVector result = NumericVector(chance.size());
  for (int i = 0; i < chance.size(); i++) {
    double individualChance = chance[i];
    result[i] = binomSniffles(x, size, individualChance);
  }
  return result;
}
