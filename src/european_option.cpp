#include "european_option.hpp"

double european_option::stock_price_single(){
    return stock_price_single(T,t);
}

double european_option::stock_price_single(double time_final, double time_inital) const{
    return St * exp((r - 0.5*sigma*sigma) * (time_final - time_inital) + sigma * sqrt(time_final - time_inital) * get_random());
}

double european_option::stock_price_single(double Z){
    return St * exp((r - 0.5*sigma*sigma) * (T - t) + sigma * sqrt(T - t) * Z);
}

double european_option::call_payoff_single(){
    double current_value = stock_price_single();
    current_value > K ? current_value = discount*(current_value - K) : current_value = 0;
    return current_value;
}

double european_option::call_payoff_single(double Z){
    double current_value = stock_price_single(Z);
    current_value > K ? current_value = discount*(current_value - K) : current_value = 0;
    return current_value;
}

double european_option::d1_calculate() {
    return (log(St/K)+(r+0.5*sigma*sigma)*tau)/sigma*sqrt(tau);
}

double european_option::d2_calculate() {
    return (log(St/K)+(r-0.5*sigma*sigma)*(tau))/sigma*sqrt(tau);
}




// MAIN FUNCTIONS

double european_option::payoff_theoretic() {
    return St*normal_cdf(d1) - K*discount*normal_cdf(d2);
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
    double result = 0.0;
    for (long int i = 0; i < number_iterations; i++){
        (stock_price_at_maturity[i] > K) ? result += discount * (stock_price_at_maturity[i] - K) * random_sample[i] : result += 0;
    }
    return result / (St * sigma * sqrt(tau) * number_iterations);
}

double european_option::delta_pw(){
    double result = 0.0;
    for (long int i = 0; i < number_iterations; i++){
        (stock_price_at_maturity[i] > K) ? result += discount * stock_price_at_maturity[i] /St : result += 0;
    }
    return result / number_iterations;
}


double european_option::gamma_lrpw(){
    double result = 0.0;
    for (long int i = 0; i < number_iterations; i++){
        (stock_price_at_maturity[i] > K) ? result +=  discount * K * random_sample[i] : result += 0;
    }
    return result / (number_iterations * St * St * sigma * sqrt(tau));
}
double european_option::gamma_pwlr(){
    double result = 0.0;
    for (long int i = 0; i < number_iterations; i++){
        (stock_price_at_maturity[i] > K) ? result +=  discount*(stock_price_at_maturity[i] / (St*St)) * (random_sample[i] / (sigma * sqrt(tau)) - 1) : result += 0;
    }
    return result / number_iterations;
}
double european_option::gamma_lrlr(){
    double result = 0.0;
    for (long int i = 0; i < number_iterations; i++){
        (stock_price_at_maturity[i] > K) ? result +=  discount*(stock_price_at_maturity[i] - K) * (( pow(random_sample[i], 2) - 1) / (St*St*sigma*sigma*(tau)) - random_sample[i] / (St*St*sigma*sqrt(tau))) : result += 0;
    }
    return result / number_iterations;
}


double european_option::vega_lr(){
    double result = 0.0;
    for (long int i = 0; i < number_iterations; i++){
        (stock_price_at_maturity[i] > K) ? result += discount * (stock_price_at_maturity[i] - K) * ((pow(random_sample[i], 2) - 1) / sigma - random_sample[i] * sqrt(tau)) : result += 0;
    }
    return result / number_iterations;
}
double european_option::vega_pw(){
    double result = 0.0;
    for (long int i = 0; i < number_iterations; i++){
        (stock_price_at_maturity[i] > K) ? result += discount * stock_price_at_maturity[i] * (-sigma*tau + sqrt(tau)*random_sample[i]) : result += 0;
    }
    return result / number_iterations;
}

// --------------
// PUBLIC
// --------------

// ACCESS FUNCTIONS

double european_option::price(){
    double result = 0.0;
    for (int i=0; i< number_iterations; i++){
        result += call_payoff_single(random_sample[i]);
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

double european_option::vega(std::string method){
    if (method=="pw") {return vega_pw();}
    if (method=="th") {return vega_theoretic();}
    return vega_lr();
}
double european_option::vega(){return vega("lr");}


// SERVICE FUNCTIONS
