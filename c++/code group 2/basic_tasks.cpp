#include <iostream>
#include <ctime> // For the random seed
#include <math.h> // For math calculations
#include <vector>

using namespace std;

// Function used to compute X->U([0,1]) sample
double uniform_random_sample(){
    int number = rand() % RAND_MAX;
    return double(number)/RAND_MAX;
}

// Function used to compute X->N(0,1) sample using Box-Muller Method
double average(vector<double> samples){
    double avg = 0;
    int i = samples.size();
    for (int j=0; j<i; j++){
        avg += samples[j]/i;
    }
    return avg;
}
double standard_deviation(vector<double> samples){
    double std = 0;
    int i = samples.size();
    double avg = average(samples);
    for (int j=0; j<i; j++){
        std += (samples[j]-avg)*(samples[j]-avg)/i;
    }
    std = sqrt(std);
    return std;
}
vector<double> normal_random_sample(int number_of_samples){
    number_of_samples>=1 ? true : number_of_samples=1;
    int iterations = 12;
    vector<double> samples;
    for(int k=0; k<number_of_samples; k++){
        double sample = 0;
        for(int i=0; i<iterations; i++){
            sample += uniform_random_sample();
        }
        samples.push_back(12*(sample - iterations*0.5)/sqrt(iterations));
    }
    double avg = average(samples); double std = standard_deviation(samples);
    for(int k=0; k<number_of_samples; k++){
       samples[k] = (samples[k] - avg)/std;
    }
    return samples;
}


// BS Closed Form version
vector<double> stock_values_closed_form(double r, double sigma, double S0, int T, int n){
    vector<double> brownian_motion = normal_random_sample(n);
    vector<double> stock_evolution;
    double S_previous = S0;
    double computed_value;
    double h = double(T)/n;
    stock_evolution.push_back(S_previous);
    for (int t=1; t<=n; t++){
        computed_value = S0*exp((r - 0.5*pow(sigma, 2)*t*h + sigma*sqrt(t*h)*brownian_motion[t]));
        S_previous = computed_value;
        stock_evolution.push_back(computed_value);
    }
    return stock_evolution;
}



// Generate stock values
vector<double> stock_values_euler_method(double r, double sigma, double S0, int T, int n){
    vector<double> brownian_motion = normal_random_sample(n);
    vector<double> stock_evolution;
    double S_previous = S0;
    double computed_value;
    double h = double(T)/n;
    double h_sqrt = sqrt(h);
    stock_evolution.push_back(S_previous);
    for (int t=1; t<=n; t++){
        computed_value = S_previous + r*h*S_previous + sigma*S_previous*h_sqrt*brownian_motion[t];
        S_previous = computed_value;
        stock_evolution.push_back(computed_value);
    }
    return stock_evolution;
}


// Functions to compute the payoff of the call
double european_call_option_payoff(double strike, double stock_value){
    double result;
    stock_value - strike < 0 ? result = 0 : result = stock_value - strike;
    return result;
}
vector<double> european_call_option_payoff_path(double strike, vector<double> stock_evolution){
    vector<double> samples = stock_evolution;
    for(int k=0; k<samples.size(); k++){
       samples[k] = european_call_option_payoff(strike, samples[k]);
    }
    return samples;
}

double contingent_claim_price(double strike, double r, int T, int t, double value_at_maturity){
    return exp(-r*(T-t))* european_call_option_payoff(strike, value_at_maturity);
}
vector<double> contingent_claim_price_paths(double strike, double r, int T, int t, vector<double> paths){
    vector<double> samples = paths;
    for(int k=0; k<samples.size(); k++){
       samples[k] = contingent_claim_price(strike, r, T, t, paths[k]);
    }
    return samples;
}





// Monte Carlo Valuation
double monte_carlo_valuation(int number_of_simualtions, double strike, double interest_rate, double sigma, double initial_value, int maturity, int number_of_steps){
    vector<double> stock_finals;
    vector<double> stock_path;
    for (int k=0; k<number_of_simualtions; k++){
        stock_path = stock_values_closed_form(interest_rate, sigma, initial_value, maturity, number_of_steps);
        stock_finals.push_back(stock_path[stock_path.size()-1]);
    }
    return average(contingent_claim_price_paths(strike, interest_rate, maturity, 0, stock_finals));
}
double monte_carlo_valuation(int number_of_simualtions, double strike, double interest_rate, double sigma, double initial_value, int maturity){return monte_carlo_valuation(number_of_simualtions, strike, interest_rate, sigma, initial_value, maturity, 12);}


int main(){
    // The random seed should only be called once to avoid resets
    srand(time(0));

    int K = 100;
    int S0 = 100;
    int T = 1;
    double sigma = 0.4;
    double r = 0.05;

    cout << monte_carlo_valuation(1000, K, r, sigma, S0, T) << endl;
    return 0;
}