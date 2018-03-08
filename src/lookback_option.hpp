#ifndef LOOKBACK_OPTION_H
#define LOOKBACK_OPTION_H

#include "random_normal.hpp" // For random normal functions
#include "european_option.hpp" // Include the european_option class for inheritance

using namespace std;

//double h = 1.0 / 1000.0, h_ = sqrt(h); //step size for the barrier

/*
 *
This class encapsulates the Black-Scholes European Lookback Option with Fixed Strike
*/

class lookback_option: public european_option{
    // VARIABLES
    double M_;  // max/min to this time
    bool c;    // true = call option ; false = put option

    // UTILITY FUNCTIONS
    // MAIN FUNCTIONS
public:

    // Constructor
    lookback_option( double initial_stock_price = 100.0, double strike = 100.0, double interest_rate = 0.05, double volatility = 0.4,double time_final_T = 1, double time_initial_t = 0, int number_iterations_approximation = 10000, double m = 150, bool cl = true): european_option(initial_stock_price, strike, interest_rate, volatility, time_final_T, time_initial_t, number_iterations_approximation){
        M_ = m;
        c = cl;
    }

    // ACCESS FUNCTIONS


    // Computes the barrier option price
    double value();


    /* Computes the delta according to the user input:
- 'pw' for pathwise method
- 'lr' for likelihood ratios method
'lr' is built as default
    */
    double delta_LR();
    double delta_PW();

    /* Computes the gamma (likelihood ratios method)  */
    double gamma_LR();

    /* Computes the vega (likelihood ratios method)  */
    double vega_LR();

    // SERVICE FUNCTIONS
};


// ??? cdf of the maximum of Brownian motion with drift

double P_M(const double &m, const double &drift, const double &stime);


// ???? pdf of the maximum of Brownian motion with drift

double rho_M(const double &m, const double &drift, const double &stime);


// ??? pdf of the maximum of Brownian motion with drift  partial differentiation wrt drift

double rho_pdrift(const double &m, const double &drift, const double &stime);


// ????  pdf of the maximum of Brownian motion with drift  partial differentiation wrt m

double rho_prime(const double &m, const double &drift, const double &stime);


// ???? newton raphson

double newtonmax(double r, double drift, double stime);



#endif //  LOOKBACK_OPTION_H
