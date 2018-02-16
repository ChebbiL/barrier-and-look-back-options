#include<iostream>
#include"Black_Scholes.h"
#define Pi 3.141592653
#define RAN_MAX 2147483647
#define luck_number 16807
using namespace std;

unsigned long int ran(){
	static unsigned long long int n = rand();
	unsigned long int m = RAN_MAX;
	n = (n*luck_number) % m;
	return n;
}
long double rann(){
	unsigned long int z = ran(), m = RAN_MAX ;
	long double x = (long double)(z ) / (long double)(m-1);
	return x;
}
double STT(double St, double K, double r, double sigma, double T, double t, double(*fun)()){
	return St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*fun());
}
double BSCall(double St, double K, double r, double sigma, double T, double t, double(*fun)()){
	double ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*fun());
	if (ST > K){
		return exp(-r*(T - t))*(ST - K);
	}
	return 0;
}
double optional_price(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)()){
	double res = 0.0;
	/*for (long int i = 0; i < M; i++){
	res = (i / (i + 1))*res + BSCall(100, 100, 0.05, 0.4, 1, 0) / (i + 1);
	}*/
	for (long int i = 0; i < M; i++){
		res = res + BSCall(100, 100, 0.05, 0.4, 1, 0,fun);
	}
	res /= M;
	return res;
}
double deltaPW(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)()){
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
double deltaLR(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)()){
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
double gammaLRPW(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)()){
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
double gammaPWLR(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)()){
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
double gammaLRLR(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)()){
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
double vegaLR(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)()){
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
double vegaPW(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)()){
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

double snd(){
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
double sndb(){
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
double snd_1(){
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
double sndb_1(){
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