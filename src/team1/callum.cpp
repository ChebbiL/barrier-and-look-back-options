using namespace std;
#include <random>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>

random_device generator;
int M = 1000; //number of paths to average over
double h = 1.0 / 1000.0, h_ = sqrt(h); //step size for the barrier

class call //basic call class, stores initial price, strike, interest rate, volatility, time to expiry
{
protected:
	double S_0, K, r, sigma, T;
public:
	call(double price, double exercise, double interest, double vol, double expiry) : S_0(price), K(exercise), r(interest), sigma(vol), T(expiry) {}
	//~call();
	/*all of the following generate the result for a single path, then use average to work out the value*/
	virtual double value();
	double delta_LR();
	double delta_PW();
	double gamma_PWLR();
	double gamma_LR();
	double vega_LR();
	virtual double valueA();
	double delta_LRA();
	double delta_PWA();
	double gamma_PWLRA();
	double gamma_LRA();
	double vega_LRA();
	//double average(double (call::*f)());
};

class barrier : public call //subclass of call for barrier options
{
protected:
	double B;
	bool up, out; //records whether barrier is up and out, down and in, etc
public:
	barrier(double price, double exercise, double interest, double vol, double expiry, double exit, bool above = true, bool outside = true) : call(price, exercise, interest, vol, expiry) { B = exit; up = above; out = outside; } //constructor for barrier. default is up and out
	double value();
	double valueA();
};

typedef double (call::*path)();

/*class dbarrier : public barrier
{
protected:
double B2;
bool down, out2;
public:
barrier(double price, double exercise, double interest, double vol, double expiry, double exit1, double exit2, bool above = true, bool outside = true)
}*/

double call::value()
{
	double mean = (r - sigma * sigma / 2.0) * T, std = sigma * sqrt(T);
	lognormal_distribution<double> L(mean, std);
	double ex = max(L(generator)*S_0 - K, 0.0);
	return ex * exp(-r * T);
}

double call::valueA()
{
	double mean = (r - sigma * sigma / 2.0) * T, std = sigma * sqrt(T);
	normal_distribution<double> L;
	double a = L(generator);
	double ex = (max(exp(a*std + mean)*S_0 - K, 0.0) + max(exp(-a * std + mean)*S_0 - K, 0.0)) / 2.0;
	return ex * exp(-r * T);
}

double call::delta_LR()
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	double instance = max(S_T - K, 0.0);
	instance *= exp(-r * T);
	instance *= Z / (S_0 * sigma * sqrt(T));
	return instance;
}

double call::delta_LRA()
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	double instance = max(S_T - K, 0.0);
	instance *= exp(-r * T);
	instance *= Z / (S_0 * sigma * sqrt(T));
	S_T = exp((-Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	double instance1 = max(S_T - K, 0.0);
	instance1 *= exp(-r * T);
	instance1 *= -Z / (S_0 * sigma * sqrt(T));
	return (instance + instance1) / 2.0;
}

double call::delta_PW() //Pathwise method of simulating delta for call - better, but can't be applied to Barrier options
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	if (S_T > K)
	{
		return (exp(-r * T)* S_T / S_0);
	}
	else
		return 0;
}

double call::delta_PWA() //antithetic above
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	double run;
	if (S_T > K)
	{
		run = (S_T / S_0);
	}
	else
		run = 0;
	S_T = exp((-Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	if (S_T > K)
	{
		run += (S_T / S_0);
	}
	return exp(-r * T)* run / 2.0;
}

double call::gamma_PWLR() //PW-LR method of calculating gamma for call
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	if (S_T > K)
	{
		return (exp(-r * T)* S_T / (S_0*S_0) * (Z / (sigma*sqrt(T)) - 1));
	}
	else
		return 0;
}

double call::gamma_PWLRA() //antithetic above
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	double run;
	if (S_T > K)
	{
		run = (S_T / (S_0*S_0) * (Z / (sigma*sqrt(T)) - 1));
	}
	else
		run = 0;
	S_T = exp((-Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	if (S_T > K)
	{
		run += (S_T / (S_0*S_0) * (-Z / (sigma*sqrt(T)) - 1));
	}
	return exp(-r * T)* run / 2.0;
}

double call::gamma_LR() //LR method of calculating gamma
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	return -(exp(-r * T)* max(S_T - K, 0.0) * ((Z*Z - 1) * sigma * sqrt(T) - Z) / (S_0 * S_0 *sigma *sigma *T));
}

double call::gamma_LRA() //antithetic above
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	double run = -(exp(-r * T)* max(S_T - K, 0.0) * ((Z*Z - 1) * sigma * sqrt(T) - Z) / (S_0 * S_0 *sigma *sigma *T));
	Z = -Z;
	S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	run += -(exp(-r * T)* max(S_T - K, 0.0) * ((Z*Z - 1) * sigma * sqrt(T) - Z) / (S_0 * S_0 *sigma *sigma *T));
	return run / 2.0;
}

double call::vega_LR() //likelihood ratio vega calculation
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	if (S_T > K) //this can be replaced by max(S_T - K, 0.0) being put in below but is this more efficient??
	{
		return (exp(-r * T) *(S_T - K) * ((Z*Z - 1) / sigma - sqrt(T) * Z));
	}
	else
		return 0;
}

double call::vega_LRA() //antithetic above
{
	normal_distribution<double> N;
	double Z = N(generator);
	double S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0, run;
	if (S_T > K) //this can be replaced by max(S_T - K, 0.0) being put in below but is this more efficient??
	{
		run = (exp(-r * T) *(S_T - K) * ((Z*Z - 1) / sigma - sqrt(T) * Z));
	}
	else
		run = 0;
	Z = -Z;
	S_T = exp((Z * sigma * sqrt(T) + (r - sigma * sigma / 2.0) * T))*S_0;
	if (S_T > K) //this can be replaced by max(S_T - K, 0.0) being put in below but is this more efficient??
	{
		run += (exp(-r * T) *(S_T - K) * ((Z*Z - 1) / sigma - sqrt(T) * Z));
	}
	return run / 2.0;
}



double average(path f, call *C) //averages values for monte carlo estimation for regular call options and greeks, not to be used for barriers. Note: call must be passed as pointer
{
	double running = 0.0;
	for (double i = 1.0; i <= double(M); i++)
	{
		double d = invoke(f, C);
		running = (1.0 - 1.0 / i)*running + d / i;
	}
	return running;
}

vector<double> sim_analysis(path f, call **C, double analytic) //runs an analysis on the estimation of an option, return vector of bias, mean error, variance(assumes there is low bias, really is MSE). Note, call must be passed as double pointer
{
	double bias = 0.0, mean = 0.0, var = 0.0;
	for (double i = 1.0; i <= 1000.0; i++)
	{
		double d = invoke(average, f, *C);
		bias = (1.0 - 1.0 / i)*bias + d / i;
		mean = (1.0 - 1.0 / i)*mean + abs(d) / i;
		var = ((i - 1.0)*var + (d - analytic)*(d - analytic)) / i;
	}
	vector<double> v;
	v.push_back(bias);
	v.push_back(mean);
	v.push_back(var);
	return v;
}


double barrier::value() //estimates the value of a barrier option
{
	double S_t = S_0, mean = (r - sigma * sigma / 2.0) * h, std = sigma * h_;
	lognormal_distribution<double> L(mean, std);
	bool hit = false;
	double time;
	for (time = 0; time <= T && !hit; time += h)
	{
		double e = L(generator);
		S_t *= e;
		if ((S_t > B) ^ !up)
			hit = true;
	}
	if (out)
	{
		if (hit)
			S_t = 0.0;
	}
	else
	{
		if (hit)
		{
			lognormal_distribution<double> L((r - sigma * sigma / 2.0) / (T - time), sigma*sqrt(T - time));
			double e = L(generator);
			S_t *= e;
		}
		else
			S_t = 0;
	}
	return max(S_t - K, 0.0)*exp(-r * T);
}

double barrier::valueA() //antithetic above
{
	double S_tp = S_0, S_tn = S_0, mean = (r - sigma * sigma / 2.0) * h, std = sigma * h_;
	normal_distribution<double> N(0.0, std);
	bool hitp = false, hitn = false;
	double time;
	for (time = 0; time < T && !(hitp && hitn); time += h)
	{
		double e = N(generator);
		S_tp *= exp(mean + e);
		S_tn *= exp(mean - e);
		if ((S_tp > B) ^ !up)
			hitp = true;
		if ((S_tn > B) ^ !up)
			hitn = true;
	}
	if (out)
	{
		if (hitp)
			S_tp = 0;
		if (hitn)
			S_tn = 0;
	}
	else
	{
		normal_distribution<double> N_(0.0, sigma*sqrt(T - time));
		double e = N_(generator);
		if (hitp)
			S_tp *= exp(mean + e);
		else
			S_tp = 0;
		if (hitn)
			S_tn *= exp(mean - e);
		else
			S_tn = 0;
	}

	return (max(S_tp - K, 0.0) + max(S_tn - K, 0.0))*exp(-r * T) / 2.0;
}

int main()
{
	barrier b(150, 120, 0.03, 0.4, 1, 200);
	call c(150, 120, 0.03, 0.4, 1);
	call *p = &b;
	call *q = &c;
	call **r = &q;
	path val = &call::valueA;
	cout << "barrier: " << average(val, s) << "\ncall: " << average(val, t) << endl;
	cout << "sim test: " << sim_analysis(val, r, 41.9344)[0];
	return 0;
}