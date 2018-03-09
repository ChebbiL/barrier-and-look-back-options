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


double lookback_option::delta(std::string method){  // choosing "pw" or "lr" method of calculation delta
    if (method=="pw") {return delta_pw();}
    return delta_lr();
}


double lookback_option::delta(){return delta("lr");}  // default mathod for delta


double lookback_option::gamma()
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


double lookback_option::vega()
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

///////////////////////////////////////

double d1(double moneyness, double r, double sigma, double T)
{
    double a = moneyness + (r + sigma*sigma / 2.0) *T;
    return a/(sigma * sqrt(T));
}

double lookback_option::Call_An()
{
    double d1 = d1(log(S_0/max(maxmin, K)), r, sigma);
    double d2 = d1 - sigma*sqrt(tau);
    double answer = St*normal_cdf(d1) - exp(-r*tau) *(normal_cdf(d2)*max(K,maxmin) + max(maxmin - K, 0.0)); 
    double ans2 = -pow(St/(max(maxmin, K)), -2*r/(sigma*sigma))*normal_cdf(d1 - r/sigma*sqrt(tau));
    ans2 += exp(r*tau)*normal_cdf(d1);
    ans2 *= St *sigma*sigma*0.5/r;
    ans2*= exp(-r*tau);
    return ans + ans2;
}

////////////////////////

#include <fstream>
#include <ctime>

int write_csv(string pathToFile){
  ofstream myfile;
  myfile.open(pathToFile);
  myfile << "iterations;delta lr;delta lr time;delta pw;delta pw time;gamma pwlr;gamma pwlr time;vega lr;vega lr time;\n";

  european_option s1(100, 100, 0.05, 0.4, 1, 0, 1);

  for (int i=1; i<10000000; i*=10){
    for (int j=0; j<1000; j+=1){
      s1 = lookback_option(100, 100, 0.05, 0.4, 1, 0, i);

   /*   clock_t begin = clock();
      double delta_th = s1.delta("th");
      clock_t end = clock();
      double delta_th_time = double(end - begin) / CLOCKS_PER_SEC;
     */ begin = clock();
      double delta_lr = s1.delta_lr;
      end = clock();
      double delta_lr_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double delta_pw = s1.delta_pw;
      end = clock();
      double delta_pw_time = double(end - begin) / CLOCKS_PER_SEC;


      //myfile << i << ";" << delta_th << ";" << delta_th_time;
      myfile << ";" << delta_lr << ";" << delta_lr_time;
      myfile << ";" << delta_pw << ";" << delta_pw_time;

     /* begin = clock();
      double gamma_th = s1.gamma("th");
      end = clock();
      double gamma_th_time = double(end - begin) / CLOCKS_PER_SEC;
    */begin = clock();
      double gamma_pwlr = s1.gamma();
      end = clock();
      double gamma_pwlr_time = double(end - begin) / CLOCKS_PER_SEC;
     /* begin = clock();
      double gamma_lrpw = s1.gamma("lrpw");
      end = clock();
      double gamma_lrpw_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double gamma_lrlr = s1.gamma("lrlr");
      end = clock();
      double gamma_lrlr_time = double(end - begin) / CLOCKS_PER_SEC;
    */

      //myfile << ";" << gamma_th << ";" << gamma_th_time;
      myfile << ";" << gamma_pwlr << ";" << gamma_pwlr_time;
      //myfile << ";" << gamma_lrpw << ";" << gamma_lrpw_time;
      //myfile << ";" << gamma_lrlr << ";" << gamma_lrlr_time;

      begin = clock();
      double vega_th = s1.vega();
      end = clock();
      double vega_th_time = double(end - begin) / CLOCKS_PER_SEC;
      /*begin = clock();
      double vega_lr = s1.vega("lr");
      end = clock();
      double vega_lr_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double vega_pw = s1.vega("pw");
      end = clock();
      double vega_pw_time = double(end - begin) / CLOCKS_PER_SEC;
*/
      myfile << ";" << vega_th << ";" << vega_th_time;
     // myfile << ";" << vega_lr << ";" << vega_lr_time;
      //myfile << ";" << vega_pw << ";" << vega_pw_time;


      myfile << "\n";

    }
  }


  return 0;
}

