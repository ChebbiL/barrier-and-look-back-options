#include "european_option.hpp"

double european_option::stock_price_single(){
    return St * exp((r - 0.5*pow(sigma, 2)) * (T - t) + sigma * sqrt(T - t) * get_random());
}

double european_option::stock_price_single(double Z){
    return St * exp((r - 0.5*pow(sigma, 2)) * (T - t) + sigma * sqrt(T - t) * Z);
}

double european_option::call_payoff_single(){
    double current_value = stock_price_single();
    current_value > K ? current_value = discount*(current_value - K) : current_value = 0;
    return current_value;
}

double european_option::d1_calculate() {
    return (log(St/K)+(r+0.5*pow(sigma,2))*(T-t))/sigma*sqrt(T-t);
}

double european_option::d2_calculate() {
    return (log(St/K)+(r-0.5*pow(sigma,2))*(T-t))/sigma*sqrt(T-t);
}




// MAIN FUNCTIONS

double european_option::payoff_theoretic() {
    return St*normal_cdf(d1) - K*exp(-r*(T-t))*normal_cdf(d2);
}

double european_option::delta_theoretic() {
    return normal_cdf(d1);
}

double european_option::gamma_theoretic() {
    return (exp(-0.5*d1*d1) / sqrt(2 * M_PI)) / (St*sigma*sqrt(T));
}

double european_option::vega_theoretic() {
    return St*sqrt(T)*(exp(-0.5*d1*d1) / sqrt(2 * M_PI));
}

double european_option::delta_lr(){
    double result = .0, Z, current_value;
    for (long int i = 0; i < number_iterations; i++){
        Z = get_random();
        current_value = stock_price_single(Z);
        (current_value > K) ? result += discount * (current_value - K) * Z : result += 0;
    }
    return result / (St * sigma * sqrt(tau) * number_iterations);
}

double european_option::delta_pw(){
    double result = 0.0, current_value;
    for (long int i = 0; i < number_iterations; i++){
        current_value = stock_price_single();
        (current_value > K / St) ? result += discount * current_value : result += 0;
    }
    return result / number_iterations;
}


double european_option::gamma_lrpw(){
    double result = .0, Z, current_value;
    for (long int i = 0; i < number_iterations; i++){
        Z = get_random();
        current_value = stock_price_single(Z);
        (current_value > K) ? result +=  discount * K * Z : result += 0;
    }
    return result/(number_iterations * St * St * sigma * sqrt(tau));
}
double european_option::gamma_pwlr(){
    double result = .0, Z, current_value;
    for (long int i = 0; i < number_iterations; i++){
        Z = get_random();
        current_value = stock_price_single(Z);
        (current_value > K) ? result +=  discount*(current_value / (St*St)) * (Z / (sigma * sqrt(tau)) - 1) : result += 0;
    }
    return result/number_iterations;
}
double european_option::gamma_lrlr(){
    double result = .0, Z, current_value;
    for (long int i = 0; i < number_iterations; i++){
        Z = get_random();
        current_value = stock_price_single(Z);
        (current_value > K) ? result +=  discount*(current_value - K) * ((Z*Z - 1) / (St*St*sigma*sigma*(tau)) - Z / (St*St*sigma*sqrt(tau))) : result += 0;
    }
    return result / number_iterations;
}

// --------------
// PUBLIC
// --------------

// ACCESS FUNCTIONS

double european_option::price(){
    double result = .0;
    for (int i=0; i< number_iterations; i++){
        result += call_payoff_single();
    }
    return result/number_iterations;
}

double european_option::delta(std::string method){
    if (method=="pw") {return delta_pw();}
    if (method=="th") {return delta_theoretic();}
    return delta_lr();
}
double european_option::delta(){return delta("lr");}

double european_option::gamma(std::string method){
    if (method=="lrpw") {return gamma_lrpw();}
    if (method=="lrlr") {return gamma_lrlr();}
    if (method=="th") {return gamma_theoretic();}
    return gamma_pwlr();
}
double european_option::gamma(){return gamma("pwlr");}



// SERVICE FUNCTIONS
