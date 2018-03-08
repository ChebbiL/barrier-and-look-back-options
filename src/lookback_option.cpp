#include "lookback_option.hpp"

// calculates the value of lookback option (average over simulations)

double lookback_option::value()
{
    double u, Z;
    long double sum = 0;
    double drift = (r / sigma - sigma / 2.0);

    for (int i=0; i<number_iterations; i++) {
        u = get_random();
        Z = newtonmax(u, drift, sqrt(tau));
        Z *= sigma;
        sum += exp(-r * tau) * (max(St * exp(Z) - max(M_, K), 0.0) + max(M_ - K, 0.0));
    }
    return sum / number_iterations;
}


// calculates the value of delta(PW) of lookback option (average over simulations)

double lookback_option::delta_PW()
{
    double u, Z;
    long double sum = 0;
    double drift = (r / sigma - sigma / 2.0);
    // loop: summing and dividing to get average
    for (int i=0; i<number_iterations; i++) {
        u = get_random();
        Z = newtonmax(u, drift, sqrt(tau));
        if ( exp(Z * sigma) * St > max(M_, K) )   // XNOR
            sum += exp(Z * sigma - r * tau);
    }
    return sum/number_iterations;
}


// calculates the value of delta(LR) of lookback option (average over simulations)

double lookback_option::delta_LR()
{
    long double sum = 0;
    double u,Z;
    double stime = sqrt(tau);
    double drift = (r / sigma - sigma / 2.0);

    for (int i=0; i<number_iterations; i++) {
        u = get_random();
        Z = newtonmax(u, drift, stime);

        if ( exp(Z * sigma) * St >= max(M_, K) ) {
            u = rho_prime(Z, drift, stime) / rho_M(Z, drift, stime);   ///// TERM
            u /= (St * sigma);                                        /////
            sum += exp(-r * tau) * u * (max(M_, K) - St * exp(Z * sigma));
        }
    }
    return sum / number_iterations;
}


// calculates the value of gamma(LR) of lookback option (average over simulations)

double lookback_option::gamma_LR()
{
    double u, Z;
    long double sum = 0;
    double stime = sqrt(tau);
    double drift = (r / sigma - sigma / 2.0);

    for (int i=0; i<number_iterations; i++) {
        u = get_random();
        Z = newtonmax(u, drift, stime);

        if ( exp(Z * sigma) * St > max(M_, K) ) {
            u = rho_prime(Z, drift, stime) / rho_M(Z, drift, stime);
            double answer = -u / (St * sigma);
            answer = answer * exp(Z * sigma);
            sum += (answer - exp(Z * sigma) / St) * exp(-r * tau);
        }
    }
    return sum / number_iterations;
}


// calculates the value of vega(LR) of lookback option (average over simulations)

double lookback_option::vega_LR()
{
    double u, Z;
    long double sum = 0;
    double temp;
    double stime = sqrt(tau);
    double drift = (r / sigma - sigma / 2.0);

    for (int i=0; i<number_iterations; i++) {
        u = get_random();
        Z = newtonmax(u, drift, stime);
        if (exp(Z * sigma) * St > max(M_, K)) {
            temp = rho_M(Z, drift, stime);
            u = temp + Z * rho_prime(Z, drift, stime);
            u = u / sigma;
            u += rho_pdrift(Z, drift, stime) * (0.5 + r / (sigma * sigma));
            u = -u / temp;
            sum += u * (St * exp(Z * sigma) - max(M_, K)) * exp(-r * tau);
        }
    }
    return sum / number_iterations;
}




// ??? cdf of the maximum of Brownian motion with drift

double P_M(const double &m, const double &drift, const double &stime)    // stime is sqrt(t)
{
    double k = drift * stime;
    return (normal_cdf(m / stime - k) - exp(2 * drift*m)*normal_cdf(-m / stime - k));
}


// ???? pdf of the maximum of Brownian motion with drift

double rho_M(const double &m, const double &drift, const double &stime)
{
    double k = drift * stime;
    double p = (2 / stime)*normal_pdf(m / stime - k);
    return p - 2 * drift * exp(2 * drift*m)*normal_cdf(-m / stime - k);
}


// ??? pdf of the maximum of Brownian motion with drift  partial differentiation wrt drift

double rho_pdrift(const double &m, const double &drift, const double &stime)
{
    double k = drift * stime;
    double x = -m / stime - k;
    double answer = -(1 + 2 * drift * m)* normal_cdf(x) / drift;   //  / drift added
    answer += normal_pdf(x)*stime;
    answer *= drift *exp(2 * drift*m);   // 2 * drift *exp(2 * drift*m)
    answer += (m / stime - k)*normal_pdf(m / stime - k);
    return answer * 2;
}


// ????  pdf of the maximum of Brownian motion with drift  partial differentiation wrt m

double rho_prime(const double &m, const double &drift, const double &stime)
{
    double k = drift * stime;
    double x = -m / stime - k;
    double answer = 2 * drift*exp(2 * drift*m)*(normal_pdf(x) / stime - 2 * drift*normal_cdf(x));
    answer -= (2 / (stime*stime))*(m/stime - k)*normal_pdf(m / stime - k);
    return answer;
}


// ???? newton raphson

double newtonmax(double r, double drift, double stime)
{
    double x = 0.0000001, c = P_M(x, drift, stime), k;
    do
    {
        k = x;
        x = x - (c - r) / rho_M(x, drift, stime);
        c = P_M(x, drift, stime);
    } while (abs(c - r) > 0.0001);
    return x;
}




