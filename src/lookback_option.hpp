#ifndef LOOKBACK_OPTION_H
#define LOOKBACK_OPTION_H

#include "european_option.hpp" // Include the european_option class for inheritance
using namespace std;

/*
 *
This class encapsulates the Black-Scholes European Lookback Call Option with Fixed Strike
*/

class lookback_option: public european_option{
    // VARIABLES
protected:
    double maxmin;   // max/min observed up to present time (t)


    // UTILITY FUNCTIONS

    // Computes the cdf of the maximum of Brownian motion with drift
    double P_M(const double &m, const double &drift, const double &stime);
    //  Computes the pdf of the maximum of Brownian motion with drift
    double rho_M(const double &m, const double &drift, const double &stime);
    // Computes the pdf of the maximum of Brownian motion with drift, partially differentiated wrt drift
    double rho_pdrift(const double &m, const double &drift, const double &stime);
    // Computes the pdf of the maximum of Brownian motion with drift, partially differentiated wrt m
    double rho_prime(const double &m, const double &drift, const double &stime);
    // Newton–Raphson method
    double newtonmax(double r, double drift, double stime);


    // MAIN FUNCTIONS

    // Computes delta using likelihood ratio method
    double delta_lr();
    // Computes delta using pathwise derivatives estimates method
    double delta_pw();

public:
    // Constructor
    lookback_option(double min_or_max_observed = 150, double initial_stock_price = 100.0, double strike = 100.0, double interest_rate = 0.05,
                     double volatility = 0.4,double time_final_T = 1, double time_initial_t = 0,
                     int number_iterations_approximation = 10000): european_option(initial_stock_price, strike, interest_rate, volatility, time_final_T, time_initial_t, number_iterations_approximation){
        maxmin = min_or_max_observed; }


    // ACCESS FUNCTIONS

    // Computes price of European Lookback Call Option with Fixed Strike
    double price();
    /* Computes the delta according to the user input:
    - 'pw' for pathwise derivatives estimates method
    - 'lr' for likelihood ratios method
    'lr' is built as default
    */
    double delta(std::string method);
    double delta();
    // Computes the gamma (likelihood ratios method)
    double gamma();
    // Computes the vega (likelihood ratios method)
    double vega();


    // SERVICE FUNCTIONS

};



#endif
