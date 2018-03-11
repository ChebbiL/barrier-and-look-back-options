#include "lookback_option.hpp"

// UTILITY FUNCTIONS


double lookback_option::P_M(const double &m, const double &drift, const double &stime) // 'stime' is sqrt(time_to_maturity)
{
    double k = drift * stime;
    return (normal_cdf(m / stime - k) - exp(2 * drift*m)*normal_cdf(-m / stime - k));
}


double lookback_option::rho_M(const double &m, const double &drift, const double &stime)
{
    double k = drift * stime;
    double p = (2 / stime)*normal_pdf(m / stime - k);
    return p - 2 * drift * exp(2 * drift*m)*normal_cdf(-m / stime - k);
}


double lookback_option::rho_pdrift(const double &m, const double &drift, const double &stime)
{
    double k = drift * stime;
    double x = -m / stime - k;
    double answer = -(1 + 2 * drift * m)* normal_cdf(x) / drift;
    answer += normal_pdf(x)*stime;
    answer *= drift *exp(2 * drift*m);
    answer += (m / stime - k)*normal_pdf(m / stime - k);
    return answer * 2;
}


double lookback_option::rho_prime(const double &m, const double &drift, const double &stime)
{
    double k = drift * stime;
    double x = -m / stime - k;
    double answer = 2 * drift*exp(2 * drift*m)*(normal_pdf(x) / stime - 2 * drift*normal_cdf(x));
    answer -= (2 / (stime*stime))*(m/stime - k)*normal_pdf(m / stime - k);
    return answer;
}


double lookback_option::newtonmax(double r, double drift, double stime) // 'stime' is sqrt(time_to_maturity)
{
    double x = 0.0001;            // initial point, close to 0 s.t. method is always convergent
    double c = P_M(x, drift, stime);
    do
    {
        x = x - (c - r) / rho_M(x, drift, stime);
        c = P_M(x, drift, stime);
    } while (abs(c - r) > 0.0001);
    return x;
}

// MAIN FUNCTIONS


double lookback_option::delta_lr()
{
    long double sum = 0;
    double u, Z, stime = sqrt(tau);;
    double drift = (r / sigma - sigma / 2.0);
    for (int i=0; i<number_iterations; i++) {  // loop: summing all values and dividing to get the average
        u = get_Urandom();
        Z = newtonmax(u, drift, stime);
        if ( exp(Z * sigma) * St >= max(maxmin, K) ) {
            u = rho_prime(Z, drift, stime) / rho_M(Z, drift, stime);
            u /= (St * sigma);
            sum += exp(-r * tau) * u * (max(maxmin, K) - St * exp(Z * sigma));
        }
    }
    return sum / number_iterations;
}


double lookback_option::delta_pw()
{
    double u, Z;
    long double sum = 0;
    double drift = (r / sigma - sigma / 2.0);
    for (int i=0; i<number_iterations; i++) {   // loop: summing all values and dividing to get the average
        u = get_Urandom();
        Z = newtonmax(u, drift, sqrt(tau));
        if ( exp(Z * sigma) * St > max(maxmin, K) )
            sum += exp(Z * sigma - r * tau);
    }
    return sum/number_iterations;
}

double lookback_option::delta_theoretic() const {
    if (St <= K) {
        double d1 = (log(St / K) + (r + sigma * sigma / 2) * T) / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);
        double d3 = d1 - 2 * r * sqrt(T) / sigma;
        double u = 1 / (sigma * sqrt(T));
        double u1 = exp(-r * T) * K / St;
        double u2 = sigma * sigma * 0.5 / r;
        double u3 = exp(-r * T) * pow(St / K, -2 * r / (sigma * sigma));
        return normal_cdf(d1) + u * normal_pdf(d1) - u1 * u * normal_pdf(d2) + u2 * (normal_cdf(d1) + u * normal_pdf(d1))
               - u3 * normal_cdf(d3) * (u2 - 1) - u3 * u2 * u * normal_pdf(d3);
    }
    else {
        double d1 = (r + sigma * sigma / 2) * T / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);
        double d3 = d1 - 2 * r * sqrt(T) / sigma;
        double u = 1 / (sigma * sqrt(T));
        double u1 = exp(-r * T);
        double u2 = sigma * sigma * 0.5 / r;
        return normal_cdf(d1) + u * normal_pdf(d1) - u1 * u * normal_pdf(d2) + u2 * (normal_cdf(d1) + u * normal_pdf(d1))
               - u1 * normal_cdf(d3) * (u2 - 1) - u1 * u2 * u * normal_pdf(d3);
    }
}


double lookback_option::gamma_theoretic() const {
    if (St <= K) {
        double d1 = (log(St / K) + (r + sigma * sigma / 2) * T) / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);
        double d3 = d1 - 2 * r * sqrt(T) / sigma;
        double u = 1 / (sigma * sqrt(T));
        double u1 = exp(-r * T) * K / St;
        double u2 = sigma * sigma * 0.5 / r;
        double u3 = exp(-r * T) * pow(St / K, -2 * r / (sigma * sigma));
        double res1 = (normal_pdf(d1)*u / St)*(1 + u2);
        double res2 = -res1*d1*u;
        double res3 = (normal_pdf(d3)*u / St)*u3*(2 - u2 + d3*u*u2);
        double res4 = (normal_cdf(d3)*u3 / St)*(1 - (1 / u2));
        double res5 = (u1*u / St)*normal_pdf(d2)*(1 + d2*u);
        return res1 + res2 + res3 + res4 + res5;
    }
    else {
        return 0;
    }
}


double lookback_option::gamma_lrlr()
{
    double u, Z, stime = sqrt(tau);
    long double sum = 0;
    double drift = (r / sigma - sigma / 2.0);
    for (int i=0; i<number_iterations; i++) {
        u = get_Urandom();
        Z = newtonmax(u, drift, stime);
        if ( exp(Z * sigma) * St > max(maxmin, K) ) {
            u = rho_prime(Z, drift, stime) / rho_M(Z, drift, stime);
            double answer = -u / (St * sigma);
            answer = answer * exp(Z * sigma);
            sum += (answer - exp(Z * sigma) / St) * exp(-r * tau);
        }
    }
    return sum / number_iterations;
}


double lookback_option::gamma_pwlr() {
    double res = 0.0, Mt = 0.0;
    double a = (r - pow(sigma, 2) / 2.0) / sigma;
    if (St > K)return 0;
    for (int i = 0; i < number_iterations; i++){
        double STT;
        STT = stock_price_single(St, T, t);
        Mt = St*exp((log(STT / St) + sqrt(pow(log(STT / St),2) - 2.0*sigma*sigma*(T - t)*log(get_Urandom()))) / 2.0);
        if (Mt < K) Mt=0.0;
        if (Mt > K){
            double f = +(2 / (sigma*Mt*sqrt(T - t)))*normal_pdf(d2_calculate(Mt, St) - 2.0 * a*sqrt(T - t));
            f -= 2.0 * a*pow(Mt / St, 2 * a / sigma)*normal_cdf(-d2_calculate(Mt, St)) / sigma / Mt;
            double df = +2.0*normal_pdf(d2_calculate(Mt, St) - 2.0*a*sqrt(T - t))*(d2_calculate(Mt, St) - 2.0*a*(T - t) - a*sqrt(T - t)) / (sigma*sigma*St*Mt*(T - t));
            df += (4.0*a*a / (sigma*sigma*Mt*St))*pow(Mt / St, 2.0*a / sigma)*normal_cdf(-d2_calculate(Mt, St));
            res += exp(-r*(T - t))*(Mt/ St)*(df/f-1/St);
        }
    }
    return res / number_iterations;


}




double lookback_option::vega_theoretic() const {
    if (St <= K) {
        double d1 = (log(St / K) + (r + sigma * sigma / 2) * T) / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);
        double d3 = d1 - 2 * r * sqrt(T) / sigma;
        double d1prime = (T*(sigma*sigma - 2 * r) - 2 * log(St / K)) / (2 * sqrt(T)*sigma*sigma);
        double d2prime = d1prime - sqrt(T);
        double d3prime = d1prime + 2 * r*sqrt(T) / (sigma*sigma);
        double p = pow(St / K, -2 * r / (sigma * sigma));
        double pprime = p*log(St / K) * 4 * r / (sigma*sigma*sigma);
        double u = sigma / r;
        double u1 = sigma * sigma * 0.5 / r;
        double u2 = exp(-r * T);
        double c1 = normal_pdf(d1)*St*d1prime*(1 + u1);
        double c2 = normal_cdf(d1)*St*u;
        double c3 = -normal_pdf(d2)*d2prime*K*u2;
        double c4 = -normal_pdf(d3)*d3prime*p*St*u1*u2;
        double c5 = -normal_cdf(d3)*St*u2*(p*u + pprime*u1);
        return (c1 + c2 + c3 + c4 + c5);
    }
    else {
        double d1 = (r + sigma * sigma / 2) * T / (sigma * sqrt(T));
        double d2 = d1 - sigma * sqrt(T);
        double d3 = d1 - 2 * r * sqrt(T) / sigma;
        double d1prime = (T*(sigma*sigma - 2 * r)) / (2 * sqrt(T)*sigma*sigma);
        double d2prime = d1prime - sqrt(T);
        double d3prime = d1prime + 2 * r*sqrt(T) / (sigma*sigma);
        double u = sigma / r;
        double u1 = sigma * sigma * 0.5 / r;
        double u2 = exp(-r * T);
        double c1 = normal_pdf(d1)*St*d1prime*(1 + u1);
        double c2 = normal_cdf(d1)*St*u;
        double c3 = -normal_pdf(d2)*d2prime*St*u2;
        double c4 = -normal_pdf(d3)*d3prime*St*u1*u2;
        double c5 = -normal_cdf(d3)*St*u2*u;
        return (c1 + c2 + c3 + c4 + c5);
    }
}


double lookback_option::vega_lr()
{
    double u, Z, stime = sqrt(tau);
    long double sum = 0;
    double temp;
    double drift = (r / sigma - sigma / 2.0);
    for (int i=0; i<number_iterations; i++) {
        u = get_Urandom();
        Z = newtonmax(u, drift, stime);
        if (exp(Z * sigma) * St > max(maxmin, K)) {
            temp = rho_M(Z, drift, stime);
            u = temp + Z * rho_prime(Z, drift, stime);
            u = u / sigma;
            u += rho_pdrift(Z, drift, stime) * (0.5 + r / (sigma * sigma));
            u = -u / temp;
            sum += u * (St * exp(Z * sigma) - max(maxmin, K)) * exp(-r * tau);
        }
    }
    return sum / number_iterations;
}


double lookback_option::vega_pw() {
    double ST, Mt,U,res=0.0,Z;
    for (int i = 1; i < number_iterations; i++){
        ST = stock_price_single(St, T, t);
        U = get_Urandom();
        Z = (log(St / St) - (r - 0.5*pow(sigma, 2))*(T - t)) / (sigma*sqrt(T - t));
        Mt = St*exp((log(ST / St) + sqrt(pow(log(ST / St), 2) - 2.0*sigma*sigma*(T - t)*log(U))) / 2.0);
        if (Mt > K){
            double A = 2.0*((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*Z)*(-sigma*(T - t) + sqrt(T - t)*Z) - 4.0*sigma*(T - t)*log(U);
            res += Mt*(-sigma*(T - t) + sqrt(T - t)*Z + 0.5*pow(pow(log(ST / St), 2) - 2.0*sigma*sigma*(T - t)*log(U), -0.5)*A);
        }
    }
    return res / number_iterations;
}


// ACCESS FUNCTIONS


double lookback_option::price()
{
    double u, Z;
    long double sum = 0;
    double drift = (r / sigma - sigma / 2.0);
    for (int i=0; i<number_iterations; i++) {   // loop: summing all values and dividing to get the average
        u = get_Urandom();
        Z = newtonmax(u, drift, sqrt(tau));
        Z *= sigma;
        sum += exp(-r * tau) * (max(St * exp(Z) - max(maxmin, K), 0.0) + max(maxmin - K, 0.0));
    }
    return sum / number_iterations;
}


double lookback_option::payoff_theoretic() {
    if (St <= K){
        return european_option::payoff_theoretic(St, K) + (St*sigma*sigma / (2 * r))*(normal_cdf(d2_calculate(St, K) + sigma*sqrt(T - t)) - discount*pow(St / K, -2 * r / pow(sigma, 2))*normal_cdf(-d2_calculate(K, St)));
    }
    else{
        return discount*(St - K) + european_option::payoff_theoretic(St, St) + (St*sigma*sigma / 2 / r)*(normal_cdf(d2_calculate(St, St) + sigma*sqrt(T - t)) - discount*normal_cdf(-d2_calculate(St, St)));
    }
}


double lookback_option::delta(std::string method){  // choosing "pw" or "lr" method of calculation delta
    if (method=="pw") {return delta_pw();}
    if (method=="th") {return delta_theoretic();}
    return delta_theoretic();
}
double lookback_option::delta(){return delta("th");}  // default mathod for delta


double lookback_option::gamma(std::string method){  // choosing "pw" or "lr" method of calculation gamma
    if (method=="lrlr") {return gamma_lrlr();}
    if (method=="pwlr") {return gamma_pwlr();}
    return gamma_theoretic();
}
double lookback_option::gamma(){return gamma("th");}  // default mathod for delta


double lookback_option::vega(std::string method){  // choosing "pw" or "lr" method of calculation vega
    if (method=="lr") {return vega_lr();}
    if (method=="pw") {return vega_pw();}
    if (method=="th") {return vega_theoretic();}
    return vega_theoretic();
}
double lookback_option::vega(){return vega("th");}  // default mathod for delta







