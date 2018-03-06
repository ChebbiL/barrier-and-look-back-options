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
  double barrier_option_price();

  // MAIN FUNCTIONS

public:
  // Constructor
  barrier_option(double barrier_value, double initial_stock_price = 100.0, double strike = 100.0, double interest_rate = 0.05, double volatility = 0.4,double time_final_T = 1, double time_initial_t = 0, int number_iterations_approximation = 10000): european_option(initial_stock_price, strike, interest_rate, volatility, time_final_T, time_initial_t, number_iterations_approximation){
    barrier = barrier_value;
  }
  // ACCESS FUNCTIONS

  // Computes the barrier option price
  double price();


  // SERVICE FUNCTIONS

};

#endif
