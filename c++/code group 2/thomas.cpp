#include <iostream>
#include <cmath> // For math calculations
#include <random>
#include <ctime> // For the random seed

#include <vector>
#include <string>

using namespace std;

// The seed for the random package is generated each time the program launches.
default_random_engine generator (time(0));


/*
This class allows you to generate easily vectors or doubles of STANDARD random normal variables.
Available methods are:
- box-muller
- marsaglia
Default method is Marsaglia (polar).
*/
class random_normal_variable{
    string method;

    // Utility
    double mean(vector<double> samples){
        double mean;
        for(int i=0; i<samples.size();i++){
            mean += samples[i];
        }
        return mean / samples.size();
    }

    double standard_deviation(vector<double> samples){
        double sd;
        double mu = mean(samples);
        for(int i=0; i<samples.size();i++){
            sd += (samples[i] - mu)*(samples[i] - mu);
        }
        return sqrt(sd / samples.size());
    }

    // Methods

    double boxmuller_sample(){
        uniform_real_distribution<double> u1(0.0, 1.0);
        uniform_real_distribution<double> u2(0.0, 1.0);
        // We return Z1 in the course formula.
        return sqrt(-2 * log(u1(generator)))*sin(2 * M_PI * u2(generator));
    }
    vector<double> boxmuller_vector(int output){
        output > 2 ? true : output = 2;
        vector<double> samples;
        for (int i=0; i<output; i++){
            samples.push_back(boxmuller_sample());
        }
        return samples;
    }
    double marsaglia_sample(){
        double v1; double v2;
        double w = 2;
        uniform_real_distribution<double> u1(0.0, 1.0);
        uniform_real_distribution<double> u2(0.0, 1.0);
        while(w>1){
            v1 = 2 * u1(generator) - 1;
            v2 = 2 * u2(generator) - 1;
            w = v1*v1 + v2*v2;
        }
        // We return N1 in the course formula.
        return sqrt((-2*log(w))/w)*v1;
    }
    vector<double> marsaglia_vector(int output){
        output > 2 ? true : output = 2;
        vector<double> samples;
        for (int i=0; i<output; i++){
            samples.push_back(marsaglia_sample());
        }
        return samples;
    }

public:
    random_normal_variable(string rv_method){
        if (rv_method=="marsaglia"){
            method = "marsaglia";
        } else if (rv_method=="box-muller"){
            method = "box-muller";
        } else{
            method = "box-muller";
        }
    }
    random_normal_variable(){ method = "marsaglia"; }

    string get_method(){ return method; }

    vector<double> sample(int size){
        if (method=="box-muller"){
            return boxmuller_vector(size);
        } else {
            return marsaglia_vector(size);
        }
    }

    double sample(){return sample(1)[0];}
    double operator() () {return sample();}

    void test(int number_samples){
        cout << "Testing random normal" << endl;
        cout << number_samples << " samples generated." << endl;
        cout << "Generation method: " << get_method() << endl;
        vector<double> these_samples = sample(number_samples);
        cout << "mean: " << mean(these_samples) << endl;
        cout << "standard deviation: " << standard_deviation(these_samples) << endl;
    }
    void test() {test(10000);}

};


/*
This class is used to compute the class price.
There are 4 method options:
- closed-form : uses the closed-form formula
- euler : uses the Euler method for path approximation in a Monte-Carlo scheme
- milstein : uses the Milstein method for path approximation in a Monte-Carlo scheme
- direct-monte-carlo : uses a direct method in a Monte-Carlo scheme
By default, there are 1000 Monte Carlo steps, and 100 steps for the Euler and Milstein schemes.
*/
class black_scholes_call_price{

    double initial_value, volatility, interest_rate, time_t, time_T, strike;
    string method;

    // Utility

    double normal_cdf(double value){
        // erfc is a standard function used to generate normal cdf. It corresponds to the error fucntion.
        return 0.5 * erfc( - value / sqrt(2) );
    }


    // Valuation

    double scheme_direct(random_normal_variable Z){
        return initial_value * exp((interest_rate - 0.5 * volatility * volatility) * (time_T - time_t) + volatility * sqrt(time_T - time_t) * Z());
    }
    double scheme_direct(double Z){
        return initial_value * exp((interest_rate - 0.5 * volatility * volatility) * (time_T - time_t) + volatility * sqrt(time_T - time_t) * Z);
    }

    double scheme_euler(random_normal_variable Z){
        int n = 100;
        double h = (time_T - time_t)/n;
        double sqh = sqrt(h);
        double S_t = initial_value;
        for (int i=0; i<n; i++){
            S_t *= 1 + interest_rate*h + volatility * sqh * Z();
        }
        return S_t;
    }

    double scheme_milstein(random_normal_variable Z){
        int n = 100;
        double h = (time_T - time_t)/n;
        double sqh = sqrt(h);
        double S_t = initial_value;
        double rv;
        for (int i=0; i<n; i++){
            rv = Z();
            S_t *= 1 + interest_rate*h + rv * ( volatility * sqh + 0.5 * h * volatility * volatility * rv) - volatility * volatility * 0.5 * h ;
        }
        return S_t;
    }

    double approximate_stock_price(random_normal_variable Z){
        if (method=="euler"){
            return scheme_euler(Z);
        } else if (method=="milstein"){
            return scheme_milstein(Z);
        } else {
            return scheme_direct(Z);
        }
    }

    double monte_carlo(int samples){
        samples > 1 ? true : samples = 1;
        random_normal_variable Z;
        double S_T;
        double C;
        for(int i=0; i<samples; i++){
            S_T = approximate_stock_price(Z);
            S_T > strike ? C += exp(-interest_rate * (time_T - time_t)) * (S_T - strike) : C += 0;
        }
        return C/samples;
    }

    double direct_valuation(){
        double tau = time_T - time_t;
        double d1 = (log(initial_value/strike) + (interest_rate + 0.5 * volatility * volatility)*tau) / (volatility * sqrt(tau));
        double d2 = d1 - volatility * sqrt(tau);
        return initial_value * normal_cdf(d1) - strike * exp( - interest_rate * tau) * normal_cdf(d2);
    }


public:
    black_scholes_call_price(double S0, double r, double sigma, double t, double T, double K, string build_method): initial_value(S0), volatility(sigma), interest_rate(r), time_t(t), time_T(T), strike(K){
        if (build_method=="euler"){
            method = "euler";
        } else if (build_method=="milstein"){
            method = "milstein";
        } else if (build_method=="direct-monte-carlo"){
            method = "direct-monte-carlo";
        } else {
            method = "closed-form";
        }
    }

    string get_method(){ return method; }
    void set_method(string new_method){ 
        if (new_method=="euler"){
            method = "euler";
        } else if (new_method=="milstein"){
            method = "milstein";
        } else if (new_method=="direct-monte-carlo"){
            method = "direct-monte-carlo";
        } else {
            method = "closed-form";
        }
    }

    double operator() (int samples){
        if (method=="euler"){
            return monte_carlo(samples);
        } else if (method=="milstein"){
            return monte_carlo(samples);
        } else if (method=="direct-monte-carlo"){
            return monte_carlo(samples);
        }else {
            return direct_valuation();
        }
    }
    double operator() () {return operator() (1000);}


    // Greeks

    /*
    In order to generate the delta, we use in the following options:
    - pw stands for pathwise derivatives estimates
    - lr stands for likelihood ratios
    LR is built as default
    */
    double delta(string estimation_method){
        random_normal_variable Z;
        int n = 10000;

        double paths_sum = 0;
        double tau = time_T - time_t;
        double discount = exp(-interest_rate * tau);
        double current_value;

        if (estimation_method == "pw"){      
            for (int i=0; i<n; i++){
                current_value = approximate_stock_price(Z);
                current_value > strike ? paths_sum += discount * current_value : paths_sum += 0;
            }
            return paths_sum/(n*initial_value);
        } else {
                double coefficient;
            for (int i=0; i<n; i++){
                current_value = approximate_stock_price(Z);
                coefficient = log(current_value / initial_value) - (interest_rate * 0.5 * volatility * volatility) * tau;
                current_value > strike ? paths_sum += discount * (current_value - strike) * coefficient : paths_sum += 0;
            }
            return paths_sum /(n * initial_value * volatility * volatility * sqrt(tau));
        } 
    }

    // Gamma is built using PW-LR estimator
    double gamma(){
        int n = 1000;

        random_normal_variable X;

        double paths_sum = 0;
        double tau = time_T - time_t;
        double discount = exp(- interest_rate * tau);
        uniform_real_distribution<double> u1(0.0, 1.0);

        for (int i=0; i<n; i++){
            approximate_stock_price(X) > strike ? paths_sum += discount * strike * u1(generator) : paths_sum += 0;
        }
        return paths_sum / (n * initial_value * initial_value * volatility * sqrt(tau));
    }

    double vega(){
        int n = 1000;

        double X;
        random_normal_variable Z;

        double paths_sum = 0;
        double tau = time_T - time_t;
        double discount = exp(- interest_rate * tau);
        double current_value;

        
        for (int i=0; i<n; i++){
            X = Z();
            current_value = scheme_direct(X);
            current_value > strike ? paths_sum += discount * current_value * X : paths_sum += 0;
        }
        return volatility * sqrt(tau) * tau * paths_sum / n;
    }

};


int main(){
    double r = 0.05;
    double sigma = 0.4;
    double T = 1.0;
    double S0 = 100.0;
    double K = 100.0;

    random_normal_variable Z;
    Z.test();
    cout << endl;

    cout << "__________________________________" << endl;
    cout << "Computing Black-Scholes Call Price" << endl;  
    black_scholes_call_price call(S0, r, sigma, 0, T, K, "euler");
    cout << call.get_method() << ": ";
    cout << call() << " (1000) - " << call(10000) << " (10,000) - " << call(100000) << " (100,000)" << endl;
    call.set_method("milstein");
    cout << call.get_method() << ": ";
    cout << call() << " (1000) - " << call(10000) << " (10,000) - " << call(100000) << " (100,000)" << endl;
    call.set_method("direct-monte-carlo");
    cout << call.get_method() << ": ";
    cout << call() << " (1000) - " << call(10000) << " (10,000) - " << call(100000) << " (100,000)" << endl;
    call.set_method("closed-form");
    cout << call.get_method() << ": " << call() << endl;
    cout << endl;

    cout << "Delta: " << call.delta("pw") << " (PW), " << call.delta("lr") << " (LR)" << endl;
    cout << "Gamma: " << call.gamma() << " (PW-LR)" << endl;
    cout << "Vega: " << call.vega() << " (PW)" << endl;
    cout << endl;



    return 0;
}