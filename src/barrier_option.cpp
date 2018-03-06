#include "barrier_option.hpp"

double barrier_option::barrier_down(){
  double current_value = St;
  double minimum_price = St;
  double h = tau / number_iterations;
  for (int i = 0; i < number_iterations; i++){
    current_value = stock_price_single(t + (i + 1) * h, t + i * h);
    if (minimum_price > current_value) minimum_price = current_value;
    if (minimum_price < barrier) return 0;
  }
  if (current_value > K) return current_value;
  return 0.0;
}

double barrier_option::barrier_up(){
  double current_value = St;
  double maximum_price = St;
  double h = tau / number_iterations;
  for (int i = 0; i < number_iterations; i++){
    current_value = stock_price_single(t + (i + 1) * h, t + i * h);
    if (maximum_price < current_value) maximum_price = current_value;
    if (maximum_price > barrier) return 0;
  }
  if (current_value > K) return current_value;
  return 0.0;
}

double barrier_option::barrier_option_price(){
  double result = .0;
  double current_value;
  for (int i = 0; i < number_iterations; i++){
    barrier < St ? current_value = barrier_down() : current_value = barrier_down();
    if (current_value!=0) result += discount * (current_value - K);
  }
  return result / number_iterations;
}

double barrier_option::price(){
  return barrier_option_price();
}
