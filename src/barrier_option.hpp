#ifndef BARRIER_OPTION_H   // header guard
#define BARRIER_OPTION_H

#include "random_normal.hpp" // For random normal functions
#include "european_option.hpp" // Include the european_option class for inheritance

/*
This class encapsultaes the Barrier Option. It inherits the european option class.
*/
class barrier_option: public european_option{
  // VARIABLES

  double barrier;
  // The base_option is the underlying european option

  // UTILITY FUNCTIONS
  double barrier_down();
  double barrier_up();

  //double payoff_theoretic(double initial_stock_price, double strike);

  double d_calculate(double x, double y);
  double d_calculate_minus(double x, double y);

  double barrier_down_theoric();
  double barrier_up_theoric();


  double delta_lr_step(double current_value);
  double gamma_lr_step(double current_value);
  double vega_lr_step(double current_value);


  // MAIN FUNCTIONS
  double barrier_option_price();

  double delta_lr();
  double delta_theoretic();

  double gamma_lr();
  double gamma_theoretic();
  double gamma_theoretic_dnorm(double x);

  double vega_lr();
  double vega_theoretic();

public:
  // Constructor
  barrier_option(double barrier_value, double initial_stock_price = 100.0, double strike = 100.0, double interest_rate = 0.05, double volatility = 0.4,double time_final_T = 1, double time_initial_t = 0, int number_iterations_approximation = 10000): european_option(initial_stock_price, strike, interest_rate, volatility, time_final_T, time_initial_t, number_iterations_approximation){
    barrier = barrier_value;
  }
  // ACCESS FUNCTIONS

  // Computes the barrier option price
  double price();
  double price_theoretic();
  /* Computes the delta according to the user input:
  - 'th' for result using closed-form formula
  - 'lr' for likelihood ratios method
  'lr' is built as default
  */
  double delta(std::string method);
  double delta();

  /* Computes the gamma according to the user input:
  - 'th' for result using closed-form formula
  - 'lr' for likelihood ratios method
  'lr' is built as default
  */
  double gamma(std::string method);
  double gamma();

  /* Computes the vega according to the user input:
  - 'th' for result using closed-form formula
  - 'lr' for likelihood ratios method
  'lr' is built as default
  */
  double vega(std::string method);
  double vega();

  // SERVICE FUNCTIONS

};

#endif
