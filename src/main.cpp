/*
This si the main project file.
*/

#include<iostream> // For console output
#include <cmath> // For math calculations
using namespace std;



/*
Black-Scholes Call Class
This class encapsultaes the BSCall.
*/
class random_normal{
  // VARIABLES
  double mean, variance;
  // UTILITY FUNCTIONS

  // MAIN FUNCTIONS

public:
  // Constructor
  random_normal(double random_variable_mean = 0, double random_variable_variance = 1){
    mean = random_variable_mean;
    variance = random_variable_variance;
  }
  // ACCESS FUNCTIONS

  // SERVICE FUNCTIONS

};


/*
Black-Scholes Call Class
This class encapsultaes the BSCall.
*/
class Black_Scholes{
  // VARIABLES
  double St,K,r,sigma,T,t;
  // UTILITY FUNCTIONS

  // MAIN FUNCTIONS

public:
  // Constructor
  Black_Scholes(double initial_stock_price = 100.0, double strike = 100.0, double interest_rate = 0.05, double volatility = 0.4,double time_final_T = 1, double time_initial_t = 0){
		St = initial_stock_price;
		K = strike;
		r = interest_rate;
		sigma = volatility;
		T = time_final_T;
		t = time_initial_t;
  }
  // ACCESS FUNCTIONS

  // SERVICE FUNCTIONS
};


// SERVICE FUNCTIONS

// Computes the cdf of the normal distribution
double normal_cdf(double x){
  return erfc( - x / sqrt(2)) / 2;
}
// Computes the pdf of the normal distribution
double normal_pdf(double x){
  return exp( - x * x / 2) / sqrt(2 * M_PI);
}


int main(){
  
  return 0;
}
