#import "random_normal.hpp"

double fun(){
  return 1.5;
}

// SERVICE FUNCTIONS

double normal_cdf(double x){
  return erfc( - x / sqrt(2)) / 2;
}

double normal_pdf(double x){
  return exp( - x * x / 2) / sqrt(2 * M_PI);
}
