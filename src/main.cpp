/*
This is the main project file.
*/

#include<iostream> // For console output
#include <cmath> // For math calculations
#include <string> // For user instructions in greeks
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
class call{
  // VARIABLES
  double St, K, r, sigma, T, t;
  double tau, discount;
  int number_iterations;

  // UTILITY FUNCTIONS

  // Computes the stock price using the direct method
  double stock_price_single(){
    return St * exp((r - 0.5*pow(sigma, 2)) * (T - t) + sigma * sqrt(T - t) * fun());
  }
  // Computes the stock price using the direct method with the random variable as an input
  double stock_price_single(double Z){
    return St * exp((r - 0.5*pow(sigma, 2)) * (T - t) + sigma * sqrt(T - t) * Z);
  }
  // Computes the call payoff using a single stock price (random)
  double call_payoff_single(){
    double current_value = stock_price_single();
  	current_value > K ? current_value = discount*(current_value - K) : current_value = 0;
  	return current_value;
  }

  // MAIN FUNCTIONS

  // Computes delta using likelihood ratio method
  double delta_lr(){
    double result = .0, Z, current_value;
  	for (long int i = 0; i < number_iterations; i++){
  		Z = fun();
  		current_value = stock_price_single(Z);
  		(current_value > K) ? result += discount * (current_value - K) * Z : result += 0;
  	}
    return result / (St * sigma * sqrt(tau) * number_iterations);
  }
  // Computes delta using pathwise derivatives estimates method
  double delta_pw(){
    double result = 0.0, current_value;
  	for (long int i = 0; i < number_iterations; i++){
  		current_value = stock_price_single();
  		(current_value > K / St) ? result += discount * current_value : result += 0;
  	}
    return result / number_iterations;
  }

public:
  // Constructor
  call(double initial_stock_price = 100.0, double strike = 100.0, double interest_rate = 0.05, double volatility = 0.4,double time_final_T = 1, double time_initial_t = 0){
		St = initial_stock_price;
		K = strike;
		r = interest_rate;
		sigma = volatility;
		T = time_final_T;
		t = time_initial_t;
    tau = T - t;
    discount = exp(-r*(T - t));
    number_iterations = 1000;
  }

  // ACCESS FUNCTIONS

  // Computes the call price by averaging (expectation) several (random) single call payoffs
  double price(){
    double result = .0;
    for (int i=0; i++; i< number_iterations){
      result += call_payoff_single();
    }
    return result/number_iterations;
  }
  /* Computes the delta according to the user input:
  - 'pw' for pathwise derivatives estimates method
  - 'lr' for likelihood ratios method
  'lr' is regarded as most accurate, so built as default
  */
  double delta(string method){
    if (method=="pw") {return delta_pw();}
    return delta_lr();
  }

  // SERVICE FUNCTIONS

  // Computes the delta
  double delta() {return delta("lr");}

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
