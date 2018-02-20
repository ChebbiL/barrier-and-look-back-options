#include<iostream>
#include<fstream>
#include"Black_Scholes1.h"
#define Pi 3.141592653
#define RAN_MAX 2147483647
#define luck_number 16807
using namespace std;

unsigned long int Black_Scholes::ran(){
	static unsigned long long int n = rand();
	unsigned long int m = RAN_MAX;
	n = (n*luck_number) % m;
	return n;
}
long double Black_Scholes::rann(){
	unsigned long int z = ran(), m = RAN_MAX ;
	long double x = (long double)(z ) / (long double)(m-1);
	return x;
}
double Black_Scholes::STT(double St, double T, double t){
	return St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*fun());
}
double Black_Scholes::BSCall(){
	double ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*fun());
	if (ST > K){
		return exp(-r*(T - t))*(ST - K);
	}
	return 0;
}
double Black_Scholes::optional_price(){
	double res = 0.0;
	/*for (long int i = 0; i < M; i++){
	res = (i / (i + 1))*res + BSCall(100, 100, 0.05, 0.4, 1, 0) / (i + 1);
	}*/
	for (long int i = 0; i < M; i++){
		res = res + BSCall();
	}
	res /= M;
	return res;
}

double Black_Scholes::deltaPW(){
	double res = 0.0, ST;
	for (long int i = 0; i < M; i++){
		ST = exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*fun());
		if (ST > K / St){
			res += exp(-r*(T - t))*ST;
		}
	}
	res /= M;
	return res;
}
double Black_Scholes::deltaLR(){
	double res = 0.0, Z, ST;
	for (long int i = 0; i < M; i++){
		Z = fun();
		ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*Z);
		if (ST > K){
			res += exp(-r*(T - t))*(ST - K)*Z / (St*sigma*sqrt(T - t));
		}
		
	}
	res /= M;
	return res;
}
double Black_Scholes::gammaLRPW(){
	double res = 0.0, Z, ST;
	for (long int i = 0; i < M; i++){
		Z = fun();
		ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*Z);
		if (ST > K){
			res += exp(-r*(T - t))*K*Z / (St*St*sigma*sqrt(T - t));
		}
	}
	res /= M;
	return res;
}
double Black_Scholes::gammaPWLR(){
	double res = 0.0, Z, ST;
	for (long int i = 0; i < M; i++){
		Z = fun();
		ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*Z);
		if (ST > K){
			res += exp(-r*(T - t))*(ST / (St*St))*(Z / (sigma*sqrt(T - t)) - 1);
		}
	}
	res /= M;
	return res;
}
double Black_Scholes::gammaLRLR(){
	double res = 0.0, Z, ST;
	for (long int i = 0; i < M; i++){
		Z = fun();
		ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*Z);
		if (ST > K){
			res += exp(-r*(T - t))*(ST - K)*((Z*Z - 1) / (St*St*sigma*sigma*(T - t)) - Z / (St*St*sigma*sqrt(T - t)));
		}
	}
	res /= M;
	return res;
}
double Black_Scholes::vegaLR(){
	double res = 0.0, Z, ST;
	for (long int i = 0; i < M; i++){
		Z = fun();
		ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*Z);
		if (ST > K){

			res += exp(-r*(T - t))*(ST - K)*((Z*Z - 1) / sigma - Z*sqrt(T - t));
		}
	}
	res /= M; 
	return res;
}
double Black_Scholes::vegaPW(){
	double res = 0.0, Z, ST;
	for (long int i = 0; i < M; i++){
		Z = fun();
		ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*Z);
		if (ST > K){
			res += exp(-r*(T - t))*ST*(-sigma*(T - t) + sqrt(T - t)*Z);
		}
	}
	res /= M;
	return res;
}

double Black_Scholes::snd(){
	static double v1, v2, W;
	static bool x=false;
	if (!x){
		do{
			v1 = 2 * rann() - 1;
			v2 = 2 * rann() - 1;
			W = v1*v1 + v2*v2;
		} while (W > 1);
		v1 *= sqrt(-2 * log(W) / W);
		v2 *= sqrt(-2 * log(W) / W);
		x = false;
		return v1;
	}
	else{
		x = true;
		return v2;
	}
}
double Black_Scholes::sndb(){
	static double v1, v2,z1,z2;
	static bool x=false;
	if (!x){
		v1 = rann();
		v2 = rann();
		z1 = sqrt(-2 * log(v1))*sin(2 * Pi * v2);
		z2 = sqrt(-2 * log(v1))*cos(2 * Pi * v2);
		x = false;
		return z1;
	}
	else{
		x = true;
		return z2;
	}
}
double Black_Scholes::snd_1(){
	static double v1, v2, W;
	static bool x = true;
	if (x){
		do{
			v1 = 2 * rann() - 1;
			v2 = 2 * rann() - 1;
			W = v1*v1 + v2*v2;
		} while (W > 1);
		v1 *= sqrt(-2 * log(W) / W);
		v2 *= sqrt(-2 * log(W) / W);
		x = false;
		return v1;
	}
	else{
		x = true;
		return v2;
	}
}
double Black_Scholes::sndb_1(){
	static double v1, v2, z1, z2;
	static bool x = true;
	if (x){
		v1 = rann();
		v2 = rann();
		z1 = sqrt(-2 * log(v1))*sin(2 * Pi * v2);
		z2 = sqrt(-2 * log(v1))*cos(2 * Pi * v2);
		x = false;
		return z1;
	}
	else{
		x = true;
		return z2;
	}
}
double Black_Scholes::normalCDF(double x)
{
	return erfc(-x / std::sqrt(2)) / 2;
}
double Black_Scholes::normalpdf(double x){
	return exp(-x*x / 2) / sqrt(2 * Pi);
}

pair<double, double> Black_Scholes::Barrier_option(double B){
	double ST = St,mt=St;
	double b = (T - t) / N;
	for (int i = 0; i < N ; i++){
		ST = STT(St, t + (i + 1)*b, t + i*b);
		if (mt > ST)mt = ST;
		if (mt < B)return make_pair(0, 0);
	}
	if (ST > K)return make_pair(ST,mt);
	else return make_pair(0, 0);
}
double Black_Scholes::Barrier_option(double B1, double B2){
	double ST = St, mt = St,Mt=St;
	double b = (T - t) / N;
	for (int i = 0; i < N; i++){
		ST = STT(St, t + (i + 1)*b, t + i*b);
		if (mt > ST)mt = ST;
		if (Mt < ST)Mt = ST;
		if (mt < B1)return 0;
		if (Mt > B2)return 0;
	}
	if (ST > K)return exp(-r*(T - t))*(ST-K);
	else return 0;
}
double Black_Scholes::EBO(double B){
	double res = 0.0;
	pair<double, double>tem;
	for (int i = 0; i < M; i++){
		tem = Barrier_option(B);
		if (tem.first != 0)res += exp(-r*(T - t))*(tem.first - K);
	}
	return res / M;
}

double Black_Scholes::EBO(double B1, double B2){
	double res = 0.0;
	for (int i = 0; i < M; i++){
		res += Barrier_option(B1,B2);
	}
	return res / M;
}
double Black_Scholes::lookback_option(double L){
	double ST = St, Mt = St;
	double b = (T - t) / N;
	for (int i = 0; i < N; i++){
		ST = STT(St, t + (i + 1)*b, t + i*b);
		if (Mt < ST)Mt = ST;
	}
	if (Mt > L)return exp(-r*(T - t))*(Mt - L);
	else return 0;
}
double Black_Scholes::ELB(double L){
	double res = 0.0;
	for (int i = 0; i < M; i++){
		res += lookback_option(L);
	}
	return res / M;
}
double Black_Scholes::BOdLR(pair<double, double> &x){
	if (x.first == 0)return 0.0;
	double mu = (r - sigma*sigma / 2) / sigma;
	double f = (1 / sqrt(T - t))*(1 / sqrt(2 * Pi))*exp(-0.5*pow(((log(x.first / St) / sigma - mu*(T-t)) / sqrt(T-t)), 2)) - exp(2 * mu*(log(x.second / St) / sigma))*(1 / sqrt(2 * Pi))*exp(-0.5*pow(((log(x.first / St) / sigma - 2 * log(x.second / St) / sigma - mu*(T-t)) / sqrt(T-t)), 2));
	double df = exp(-0.5*pow((log(x.first / St) / sigma - mu*(T - t)), 2) / (T - t))*(log(x.first / St) / sigma - mu*(T - t)) / (sqrt(2 * Pi*pow((T - t), 3))*sigma*St);
	df += sqrt(2 / Pi) *mu*x.second*pow(x.second / St, (2 * mu / sigma - 1))*exp(-pow((-2 * log(x.second / St) / sigma + log(x.first / St) / sigma - mu*(T - t)), 2) / (2 * (T - t))) / (sigma*St*St*sqrt(T - t));
	df += pow(x.second / St, (2 * mu / sigma))*(-2 * log(x.second / St) / sigma + log(x.first / St) / sigma - mu*(T - t))*exp(-pow((-2 * log(x.second / St) / sigma + log(x.first / St) / sigma - mu*(T - t)), 2) / (2 * (T - t))) / (sqrt(2 * Pi*pow((T - t), 3))*sigma*St);
	double res = exp(-r*(T - t))*(x.first - K)*df / (f*sigma*sigma*x.first*x.second);
	if (res< 0){
		ofstream gout;
		gout.open("res.txt", ios::app);
		gout << res << endl;
		gout.close();
	}
	return res;
}
double Black_Scholes::BOdeltaLR(double B){
	double res = 0.0;
	for (int i = 0; i < M; i++){
		res += BOdLR(Barrier_option(B));
	}
	return res / M;
}
double Black_Scholes::BOdeltaCF(double B){
	double p = 0.5 - r / (sigma*sigma);
	double d2 = (log(St / K) + (r - sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
	double d3 = (log(St / B) - (r - sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
	return exp(-r*(T - t))*(normalCDF(-d2) + pow(St / B, 2 * p)*normalCDF(-d3));
}

void Errortest(Black_Scholes x,double(Black_Scholes::*foo)(), double y,ostream &os){
	double ermean = 0, erms = 0;
	double tem;
	for (int i = 0; i < 100; i++){
		tem = (x.*foo)() - y;
		ermean += tem;
		erms += tem*tem;
	}
	os  << (ermean /= 100 )<< " , " << (erms / 100 - ermean*ermean);
}

