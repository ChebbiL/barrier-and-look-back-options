#include<iostream>
#include<string>
#include<thread>
#include<fstream>
#include<cmath>
#include<random>
#include<mutex> 
#include"Black_Scholes1.h"
#define Pi 3.141592653
#define RAN_MAX 2147483647
#define luck_number 16807
using namespace std;
random_device rd;
mt19937 mt(rd());
uniform_real_distribution<double> ud(0, 1.0);

unsigned long int Black_Scholes::ran()const{
	static unsigned long long int n = rand();
	//unsigned long int m = RAN_MAX;
	n = (n*luck_number) % RAN_MAX;
	return n;
}
long double Black_Scholes::rann()const{
	if (method == 0){
	unsigned long int z = ran(), m = RAN_MAX ;
	long double x = (long double)(z ) / (long double)(m-1);
	return x;
	}
	else if (method == 1){
		return ud(mt);
	}
}
double Black_Scholes::BSC(double S0, double k)const{
	double d = (log(S0 / k) + (r + sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
	return normalCDF(d)*S0 - normalCDF(d - sigma*sqrt(T - t))*k*exp(-r*(T - t));
}
double Black_Scholes::BSCdelta(double S0, double k)const{
	double d = (log(S0 / k) + (r + sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
	return  normalCDF(d);
}
double Black_Scholes::BSCgamma(double S0, double k)const{
	double d = (log(S0 / k) + (r + sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
	return  normalpdf(d) / (S0*sigma*sqrt(T - t));
}
double Black_Scholes::BSCvega(double S0, double k)const{
	double d = (log(S0 / k) + (r + sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
	return  S0*normalpdf(d)*sqrt(T - t);
}
double Black_Scholes::STT(double S0, double T1, double t1)const{
	return S0*exp((r - 0.5*pow(sigma, 2))*(T1 - t1) + sigma*sqrt(T1 - t1)*fun());
}
double Black_Scholes::getmax(double S0, double SE, double T, double t)const{
	return S0*exp((log(SE / S0) + sqrt(pow(log(SE / S0),2) - 2.0*sigma*sigma*(T - t)*log(rann()))) / 2.0);
}
double Black_Scholes::getmin(double S0, double SE, double T, double t)const{
	return S0*exp((log(SE / S0) - sqrt(pow(log(SE / S0), 2) - 2.0*sigma*sigma*(T - t)*log(rann()))) / 2.0);
}

Black_Scholes::~Black_Scholes(){
	
}
double Black_Scholes::OP()const{
	double res = 0.0,ST;
	for (long int i = 0; i < M; i++){
		ST = St*exp((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*fun());
		if (ST > K){
			res+= exp(-r*(T - t))*(ST - K);
		}
	}
	res /= M;
	return res;
}

double Black_Scholes::deltaPW()const{
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
double Black_Scholes::deltaLR()const{
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
double Black_Scholes::gammaLRPW()const{
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
double Black_Scholes::gammaPWLR()const{
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
double Black_Scholes::gammaLR()const{
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
double Black_Scholes::vegaLR()const{
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
double Black_Scholes::vegaPW()const{
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

double Black_Scholes::snd()const{
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
double Black_Scholes::sndb()const{
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
double Black_Scholes::snd_1()const{
	//mtx.lock();
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
		//mtx.unlock();
		return v1;
	}
	else{
		x = true;
		//mtx.unlock();
		return v2;
	}
}
double Black_Scholes::sndb_1()const{
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
double Black_Scholes::normalCDF(double x)const
{
	return erfc(-x / sqrt(2.0)) / 2.0;
}
double Black_Scholes::normalpdf(double x)const{
	return exp(-x*x / 2.0) / sqrt(2.0 * Pi);
}
double Black_Scholes::dbs(double x, double y)const{
	return (log(x / y) + (r - sigma*sigma / 2.0)*(T - t)) / (sigma*sqrt(T - t));
}







double Barrier_option::Barrier_option_down()const{
	double ST = St,mt=St;
	//double b = (T - t) / N;
	//for (int i = 0; i < N ; i++){
	//	ST = STT(ST, t + (i + 1)*b, t + i*b);
	//	if (mt > ST)mt = ST;
	//	if (mt < B)return 0;
	//}
	//if (ST > K)return ST;
	//else return 0;
	ST = STT(St, T, t);
	mt = getmin(St, ST, T, t);
	if (mt>B && ST>K)return ST;
	else return 0.0;
}
double Barrier_option::Barrier_option_up()const{
	double ST = St, Mt = St;
	/*double b = (T - t) / N;
	for (int i = 0; i < N; i++){
		ST = STT(ST, t + (i + 1)*b, t + i*b);
		if (Mt < ST)Mt = ST;
		if (Mt > B)return 0.0;
	}
	if (ST > K)return ST;	
	else return 0.0;*/
	ST = STT(St, T, t);
	Mt = getmax(St, ST, T, t);
	if (Mt<B && ST>K)return ST;
	else return 0.0;
}
double Barrier_option::OP()const{
	double res = 0.0;
	double tem;
	for (int i = 0; i < M; i++){
		if(B<St)tem = Barrier_option_down();
		else tem = Barrier_option_up();
		if (tem!= 0)res += exp(-r*(T - t))*(tem - K);
	}
	return res / M;
}



double Barrier_option::BOdLR(double x)const{
	if (x == 0 )return 0.0;
	double mu = (r - sigma*sigma / 2) / sigma;
	double f = (normalpdf((log(x / St) - sigma*mu*(T - t)) / (sigma*sqrt(T - t))) - exp(2 * mu*log(B / St) / sigma)*normalpdf((log(x / St) - 2 * log(B / St) - sigma*mu*(T - t)) / (sigma*sqrt(T - t)))) / sqrt(T - t);
	double df = exp(-0.5*pow((log(x / St) / sigma - mu*(T - t)), 2) / (T - t))*(log(x / St) / sigma - mu*(T - t)) / (sqrt(2 * Pi*pow((T - t), 3))*sigma*St);
	df += sqrt(2 / Pi) *mu*B*pow(B / St, (2 * mu / sigma - 1))*exp(-pow((-2 * log(B / St) / sigma + log(x / St) / sigma - mu*(T - t)), 2) / (2 * (T - t))) / (sigma*St*St*sqrt(T - t));
	df += pow(B / St, (2 * mu / sigma))*(-2 * log(B / St) / sigma + log(x / St) / sigma - mu*(T - t))*exp(-pow((-2 * log(B / St) / sigma + log(x / St) / sigma - mu*(T - t)), 2) / (2 * (T - t))) / (sqrt(2 * Pi*pow((T - t), 3))*sigma*St);
	double res = exp(-r*(T - t))*(x - K)*df / f;
	return res;
}
double Barrier_option::deltaLR()const{
	double res = 0.0;
	if (B<St){
		for (int i = 0; i < M; i++){
			res += BOdLR(Barrier_option_down());
		}
	}
	else{
		for (int i = 0; i < M; i++){
			res += BOdLR(Barrier_option_up());
		}
	}
	return res / M;
}
double Barrier_option::BOCF_down()const{
	double mu = r - sigma*sigma / 2.0;
	
	double res;
	if (B <= K){
		res = BSC(St, K) - pow(St / B, -2.0 * mu / (sigma*sigma))*BSC(B*B / St, K);
	}
	else{
		res = BSC(St, B) - pow(St / B, -2.0 * mu / (sigma*sigma))*BSC(B*B / St, B);
		res += (B - K)*exp(-r*(T - t))*(normalCDF(dbs(St, B)) - pow(B / St, 2 * mu / pow(sigma, 2))*normalCDF(dbs(B*B / St, B)));
	}
	return res;
}
double Barrier_option::BOCF_up()const{
	if (B < K)return 0;
	double mu = r - sigma*sigma / 2.0;
	return BSC(St, K) - BSC(St, B) - (B - K)*exp(-r*(T - t))*normalCDF(dbs(St, B)) - pow(B / St, 2.0 * mu / (sigma*sigma))*(BSC(B*B / St, K) - BSC(B*B / St, B) - (B - K)*exp(-r*(T - t))*normalCDF(dbs(B, St)));
}
double Barrier_option::CF()const{
	if (B < St)return BOCF_down();
	else return BOCF_up();
}
double dnorm(double x){
	return -exp(-0.5*x*x)*x / sqrt(2.0*Pi);
}
double Barrier_option::deltaCF()const{
	double mu = r - sigma*sigma / 2.0;
	double res = 0.0;
	if (B < St){
		res = normalCDF((log(St / K) + (r + sigma*sigma / 2.0)*(T - t)) / (sigma*sqrt(T - t)));
		res -= pow((B / St), (r / (sigma*sigma) - 1))*(-B*B* normalCDF((log(B*B / (St*K)) + mu*(T - t)) / (sigma*sqrt(T - t)) + sigma*sqrt(T - t)) / (St*St) - 2.0 * mu*BSC(B*B / St, K) / (St*sigma*sigma));
		
	}
	else{
		res = BSCdelta(St, K) - BSCdelta(St, B);
		res -= (B - K)*exp(-r*(T - t))*normalpdf(dbs(St,B))/(sigma*St*sqrt(T-t));
		res += (2.0*mu*pow(B / St, 2.0*mu / pow(sigma, 2)) / (pow(sigma, 2)*St))*(BSC(B*B / St, K) - BSC(B*B / St, B) - (B - K)*exp(-r*(T - t))*normalCDF(dbs(B, St)));
		res -= pow(B / St, 2.0*mu / pow(sigma, 2))*(-B*B / (St*St)*BSCdelta(B*B / St, K) + B*B / (St*St)*BSCdelta(B*B / St, B) + (B - K)*exp(-r*(T - t))*normalpdf(dbs(B, St)) / (sigma*St*sqrt(T - t)));
	}
	return res;
}
double Barrier_option::BOgLR(double x)const{
	if (x == 0)return 0.0;
	double mu = (r - sigma*sigma / 2.0) / sigma;
	double f = (normalpdf((log(x / St) - sigma*mu*(T - t)) / (sigma*sqrt(T - t))) - exp(2 * mu*log(B / St) / sigma)*normalpdf((log(x / St) - 2 * log(B / St) - sigma*mu*(T - t)) / (sigma*sqrt(T - t)))) / sqrt(T - t);
	double df = exp(-pow((log(x / St) / sigma - mu*(T - t)), 2) / (2.0*(T - t)))*pow(log(x / St) / sigma - mu*(T-t), 2) / (sigma*sigma*St*St*pow((T-t), 2.5)*sqrt(2*Pi));
	df -= exp(-pow(log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T - t))) / (sqrt(2.0 * Pi*pow((T - t), 3))*sigma*sigma*St*St);
	df -= exp(-pow(log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T - t)))*(log(x / St) / sigma - mu*(T-t)) / (sqrt(2.0 * Pi*pow((T - t), 3))*sigma*St*St);
	df -= sqrt(2.0 / Pi)*mu*B*B*(2.0 * mu / sigma - 1.0)*pow(B / St, 2.0 * mu / sigma - 2.0)*exp(-pow(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T - t))) / (sigma*pow(St, 4)*sqrt(T-t));
	df -= sqrt(2.0 / Pi)*mu*B*pow(B / St, 2.0*mu / sigma - 1.0)*(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t))*exp(-pow(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T - t))) / (sigma*sigma*pow(St, 3)*pow(T - t, 1.5));
	df -= 2.0*sqrt(2.0 / Pi)*mu*B*pow(B / St, 2 * mu / sigma - 1.0)*exp(-pow(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T-t))) / (sigma*pow(St, 3)*sqrt(T - t));
	df -= sqrt(2.0 / Pi)*mu*B*pow(B / St, 2.0*mu / sigma - 1.0)*(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t))*exp(-pow(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T - t))) / (sigma*sigma*pow(T - t, 1.5)*pow(St, 3));
	df -= pow(B / St, 2 * mu / sigma)*pow(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t), 2)*exp(-pow(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T - t))) / (sqrt(2 * Pi*pow(T - t, 5))*sigma*sigma*St*St);
	df += pow(B / St, 2 * mu / sigma)*exp(-pow(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T - t))) / (sqrt(2 * Pi*pow(T - t, 3))*sigma*sigma*St*St);
	df -= pow(B / St, 2 * mu / sigma)*(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t))*exp(-pow(-2.0*log(B / St) / sigma + log(x / St) / sigma - mu*(T - t), 2) / (2.0*(T - t))) / (sqrt(2 * Pi*pow(T - t, 3))*sigma*St*St);
	return exp(-r*(T - t))*(x - K)*df / f;
}
double Barrier_option::gammaLR()const{
	double res = 0.0;
	if (B < St){
		for (int i = 0; i < M; i++){
			res += BOgLR(Barrier_option_down());
		}
	}
	else{
		for (int i = 0; i < M; i++){
			res += BOgLR(Barrier_option_up());
		}
	}
	
	return res / M;
}
double Barrier_option::gammaCF()const{
	if (B < St){
		double mu = r - sigma*sigma / 2.0;
		double d = (log(St / K) + (r + sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
		double dbs = (log(B*B / (St*K)) + mu*(T - t)) / (sigma*sqrt(T - t)) + sigma*sqrt(T - t);
		double delta = -B*B / (St*St)*normalCDF(dbs);
		double gamma = (B*B / pow(St, 3))*(2.0*normalCDF(dbs) + normalpdf(dbs) / (sigma*sqrt(T - t)));
		double res = normalpdf(d) / (St*sigma*sqrt(T - t));
		res -= pow(B / St, 2.0*mu / (sigma*sigma))*(2.0*mu*(2.0*mu + sigma*sigma)*BSC(B*B / St, K) / (St*St*pow(sigma, 4)) + gamma - 4.0*mu / (St*sigma*sigma)*delta);
		return res;
	}
	else {
		double mu = r - sigma*sigma / 2.0;
		double res = BSCgamma(St, K) - BSCgamma(St, B) - (B - K)*exp(-r*(T - t))*dnorm(dbs(St, B)) / (sigma*sigma*St*St*(T - t)) + (B - K)*exp(-r*(T - t))*normalpdf(dbs(St, B)) / (sigma*St*St*sqrt(T - t));
		res -= (4.0*mu*(mu + 0.5*sigma*sigma)*pow(B / St, 2.0*mu / sigma*sigma) / (pow(sigma, 4)*St*St))*(BSC(B*B / St, K) - BSC(B*B / St, B) - (B - K)*exp(-r*(T - t))*normalCDF(dbs(B, St)));
		res += (4.0*mu*pow(B / St, 2.0*mu / pow(sigma, 2)) / (pow(sigma, 2)*St))*(-B*B / (St*St)*BSCdelta(B*B / St, K) + B*B / (St*St)*BSCdelta(B*B / St, B) + (B - K)*exp(-r*(T - t))*normalpdf(dbs(B, St)) / (sigma*St*sqrt(T - t)));
		res -= pow(B / St, 2.0*mu / pow(sigma, 2))*(pow(B / St, 4)*BSCgamma(B*B / St, K) +2.0*B*B/pow(St,3)*BSCdelta(B*B/St,K)- pow(B / St, 4)*BSCgamma(B*B / St, B)-2.0*B*B/pow(St,3)*BSCdelta(B*B/St,B) - (B - K)*exp(-r*(T - t))*dnorm(dbs(B, St)) / (sigma*sigma*St*St*(T - t)) - (B - K)*exp(-r*(T - t))*normalpdf(dbs(B, St)) / (sigma*St*St*sqrt(T - t)));
		return res;
	}
}

double Barrier_option::BOvLR(double x)const{
	if (x == 0 )return 0.0;
	double mu = (r - sigma*sigma / 2) / sigma;
	double f = (normalpdf((log(x / St) - sigma*mu*(T - t)) / (sigma*sqrt(T - t))) - exp(2 * mu*log(B / St) / sigma)*normalpdf((log(x / St) - 2 * log(B / St) - sigma*mu*(T - t)) / (sigma*sqrt(T - t)))) / sqrt(T - t);
	f /= (sigma*x);
	double df = -exp(-pow(log(x / St) / sigma - (T - t)*(r / sigma - sigma / 2.0), 2) / (2.0*(T - t)))*((T - t)*(r / pow(sigma, 2) + 0.5) - log(x / St) / pow(sigma, 2))*(log(x / St) / sigma - (T - t)*(r / sigma - sigma / 2.0)) / (sqrt(2.0*Pi*pow(T - t, 3))*sigma*x);
	df -= exp(-pow(log(x / St) / sigma - (T - t)*(r / sigma - sigma / 2.0), 2) / (2.0*(T - t))) / (sqrt(2.0*Pi*(T - t))*pow(sigma, 2)*x);
	df += (4.0*r*log(B / St) / pow(sigma, 3) + (2.0*log(B / St) / pow(sigma, 2) + (T - t)*(r / pow(sigma, 2) + 0.5) - log(x / St) / pow(sigma, 2))*(-2.0*log(B / St) / sigma + (T - t)*(-r / sigma + sigma / 2.0) + log(x / St) / sigma) / (T - t))*exp(2.0*(r / pow(sigma, 2) - 0.5)*log(B / St) - pow(-2.0*log(B / St) / sigma + (T - t)*(-r / sigma + sigma / 2.0) + log(x / St) / sigma, 2) / (2.0*(T - t))) / (sqrt(2.0*Pi*(T - t))*sigma*x);
	df += exp(2.0*(r / pow(sigma, 2) - 0.5)*log(B / St) - pow(-2.0*log(B / St) / sigma + (T - t)*(-r / sigma + sigma / 2.0) + log(x / St) / sigma, 2) / (2.0*(T - t))) / (sqrt(2.0*Pi*(T - t))*pow(sigma, 2)*x);
	return exp(-r*(T - t))*(x - K)*df / f;
}


double Barrier_option::vegaLR()const{
	double res = 0.0;
	if (B < St){
		for (int i = 0; i < M; i++){
			res += BOvLR(Barrier_option_down());
		}
	}
	else{
		for (int i = 0; i < M; i++){
			res += BOvLR(Barrier_option_up());
		}
	}
	return res / M;
}
double Barrier_option::vegaCF()const{
	if (B < St){
		double mu = r - sigma*sigma / 2.0;
		double d = (log(St / K) + (r + sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
		double dbs = (log(B*B / (St*K)) + mu*(T - t)) / (sigma*sqrt(T - t));
		double vega = sqrt(T - t)*K*exp(-r*(T - t))*normalpdf(dbs);
		double res = St*normalpdf(d)*sqrt(T - t);
		res -= pow(B / St, 2.0*mu / (sigma*sigma))*(vega - 4.0*r*BSC(B*B / St, K)*log(B / St) / pow(sigma, 3));
		return res;
	}
	else{
		double mu = r - sigma*sigma / 2.0;
		double res = BSCvega(St, K) - BSCvega(St, B) + (B - K)*exp(-r*(T - t))*normalpdf(dbs(St, B))*((T - t)*(sigma*sigma + mu) + log(St / B)) / (sigma*sigma*sqrt(T - t));
		res += 4.0*r*log(B / St)*pow(B / St, 2 * r / pow(sigma, 2) - 1) / pow(sigma, 3)*(BSC(B*B / St, K) - BSC(B*B / St, B) - (B - K)*exp(-r*(T - t))*normalCDF(dbs(B, St)));
		res -= pow(B / St, 2.0*mu / pow(sigma, 2))*(sqrt(T - t)*exp(-r*(T - t))*(K*normalpdf(dbs(B*B / St, K)) - B*normalpdf(dbs(B*B / St, B))) + (B - K)*exp(-r*(T - t))*normalpdf(dbs(B, St))*((T - t)*(sigma*sigma + mu) + log(B / St)) / (sigma*sigma*sqrt(T - t)));
		return res;
	}
	return 0.0;
}
void Black_Scholes::input(string x){
	string tem = x,tem2;
	double num;
	int i = 0;
	while (i < tem.size()){
		
		tem2 = "";
		for (; i < tem.size() && tem[i] == 32; i++);
		for (; i < tem.size() && tem[i] != '='; i++){
			if (tem[i] != 32){
				tem2 += tem[i];
			}
		}
		for (; i < tem.size() && (tem[i]<48 || tem[i]>57); i++);
		if (i<x.size()){
			num = stod(tem.substr(i));
			for (; i < tem.size() && tem[i] != ','; i++);
			if (tem2 == "St"){ St = num; }
			else if (tem2 == "K"){ K = num; }
			else if (tem2 == "r"){ r = num; }
			else if (tem2 == "sigma"){ sigma = num; }
			else if (tem2 == "T"){ T = num; }
			else if (tem2 == "t"){ t = num; }
			else if (tem2 == "M"){ M = num; }
			i++;
		}
		
	}
}
void Barrier_option::input(string x){
	string tem = x, tem2;
	double num;
	int i = 0;
	while (i < tem.size()){
		tem2 = "";
		for (; i < tem.size() && tem[i] != '='; i++){
			if (tem[i] != 32){
				tem2 += tem[i];
			}
		}
		for (; i < tem.size() && (tem[i]<48 || tem[i]>57); i++);
		if (i<x.size()){
			num = stod(tem.substr(i));
			for (; i < tem.size() && tem[i] != ','; i++);
			if (tem2 == "St"){ St = num; }
			else if (tem2 == "K"){ K = num; }
			else if (tem2 == "r"){ r = num; }
			else if (tem2 == "sigma"){ sigma = num; }
			else if (tem2 == "T"){ T = num; }
			else if (tem2 == "t"){ t = num; }
			else if (tem2 == "B"){ B = num; }
			else if (tem2 == "M"){ M = num; }
			else if (tem2 == "N"){ N = num; }
			i++;
		}

	}
}
void Lookback_option::input(string x){
	string tem = x, tem2;
	double num;
	int i = 0;
	while (i < tem.size()){
		tem2 = "";
		for (; i < tem.size() && tem[i] != '='; i++){
			if (tem[i] != 32){
				tem2 += tem[i];
			}
		}
		for (; i < tem.size() && (tem[i]<48 || tem[i]>57); i++);
		if (i<x.size()){
			num = stod(tem.substr(i));
			for (; i < tem.size() && tem[i] != ','; i++);
			if (tem2 == "St"){ St = num; }
			else if (tem2 == "K"){ K = num; }
			else if (tem2 == "r"){ r = num; }
			else if (tem2 == "sigma"){ sigma = num; }
			else if (tem2 == "T"){ T = num; }
			else if (tem2 == "t"){ t = num; }
			else if (tem2 == "M"){ M = num; }
			else if (tem2 == "N"){ N = num; }
			i++;
		}

	}
}

double Lookback_option::lookback_option()const{
	double ST = St, Mt = St;
	double b = (T - t) / N;
	/*for (int i = 0; i < N; i++){
		ST = STT(ST, t + (i + 1)*b, t + i*b);
		if (Mt < ST)Mt = ST;
	}
	if (Mt > K)return Mt;*/
	ST = STT(St, T, t);
	Mt = getmax(St, ST, T, t);
	if (Mt > K)return Mt;
	else return 0.0;
}
double Lookback_option::OP()const{
	double res = 0.0, tem = 0.0;
	for (int i = 0; i < M; i++){
		tem = lookback_option();
		if (tem > K)res += exp(-r*(T - t))*(tem - K);
	}
	return res / M;
}
double Lookback_option::CF()const{
	if (St <= K){
		return BSC(St, K) + (St*sigma*sigma / (2 * r))*(normalCDF(dbs(St, K) + sigma*sqrt(T - t)) - exp(-r*(T - t))*pow(St / K, -2 * r / pow(sigma, 2))*normalCDF(-dbs(K, St)));
	}
	else{
		return exp(-r*(T - t))*(St - K) + BSC(St, St) + (St*sigma*sigma / 2 / r)*(normalCDF(dbs(St, St) + sigma*sqrt(T - t)) - exp(-r*(T - t))*normalCDF(-dbs(St, St)));
	}
	
}
double Lookback_option::deltaCF()const{
	if (St <= K) {
		double d1 = (log(St / K) + (r + sigma * sigma / 2) * T) / (sigma * sqrt(T));
		double d2 = d1 - sigma * sqrt(T);
		double d3 = d1 - 2 * r * sqrt(T) / sigma;
		double u = 1 / (sigma * sqrt(T));
		double u1 = exp(-r * T) * K / St;
		double u2 = sigma * sigma * 0.5 / r;
		double u3 = exp(-r * T) * pow(St / K, -2 * r / (sigma * sigma));
		return normalCDF(d1) + u * normalpdf(d1) - u1 * u * normalpdf(d2) + u2 * (normalCDF(d1) + u * normalpdf(d1))
			- u3 * normalCDF(d3) * (u2 - 1) - u3 * u2 * u * normalpdf(d3);
	}
	else {
		double d1 = (r + sigma * sigma / 2) * T / (sigma * sqrt(T));
		double d2 = d1 - sigma * sqrt(T);
		double d3 = d1 - 2 * r * sqrt(T) / sigma;
		double u = 1 / (sigma * sqrt(T));
		double u1 = exp(-r * T);
		double u2 = sigma * sigma * 0.5 / r;
		return normalCDF(d1) + u * normalpdf(d1) - u1 * u * normalpdf(d2) + u2 * (normalCDF(d1) + u * normalpdf(d1))
			- u1 * normalCDF(d3) * (u2 - 1) - u1 * u2 * u * normalpdf(d3);
	}
}
double Lookback_option::gammaCF()const{
	if (St <= K) {
		double d1 = (log(St / K) + (r + sigma * sigma / 2) * T) / (sigma * sqrt(T));
		double d2 = d1 - sigma * sqrt(T);
		double d3 = d1 - 2 * r * sqrt(T) / sigma;
		double u = 1 / (sigma * sqrt(T));
		double u1 = exp(-r * T) * K / St;
		double u2 = sigma * sigma * 0.5 / r;
		double u3 = exp(-r * T) * pow(St / K, -2 * r / (sigma * sigma));
		double res1 = (normalpdf(d1)*u / St)*(1 + u2);
		double res2 = -res1*d1*u;
		//double res3 = (normalPDF(d3)*u/S)*u3*(2-2*u2+u2*d3*u); 
		//double res4 = (normalCDF(d3)*u3/S)*(1-(1/u2)); 
		double res3 = (normalpdf(d3)*u / St)*u3*(2 - u2 + d3*u*u2);
		double res4 = (normalCDF(d3)*u3 / St)*(1 - (1 / u2));
		double res5 = (u1*u / St)*normalpdf(d2)*(1 + d2*u);
		return res1 + res2 + res3 + res4 + res5;
	}
	else {
		return 0;
	}
}
double Lookback_option::vegaCF()const{
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
		double c1 = normalpdf(d1)*St*d1prime*(1 + u1);
		double c2 = normalCDF(d1)*St*u;
		double c3 = -normalpdf(d2)*d2prime*K*u2;
		double c4 = -normalpdf(d3)*d3prime*p*St*u1*u2;
		double c5 = -normalCDF(d3)*St*u2*(p*u + pprime*u1);
		return (c1 + c2 + c3 + c4 + c5) ;
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
		double c1 = normalpdf(d1)*St*d1prime*(1 + u1);
		double c2 = normalCDF(d1)*St*u;
		double c3 = -normalpdf(d2)*d2prime*St*u2;
		double c4 = -normalpdf(d3)*d3prime*St*u1*u2;
		double c5 = -normalCDF(d3)*St*u2*u;
		return (c1 + c2 + c3 + c4 + c5) ;
	}
}
double Lookback_option::deltaPW()const{
	double res = 0.0, tem = 0.0;
	for (int i = 0; i < M; i++){
		tem = lookback_option();
		if (tem > K)res += exp(-r*(T - t))*(tem / St);
	}
	return res / M;
}
double min(double x, double y){ if (x > y)return y; return x; }
double Lookback_option::deltaLR()const{
	double res = 0.0,Mt;
	double a = (r - pow(sigma, 2) / 2.0) / sigma;
	double mu = r - pow(sigma, 2) / 2.0;
	for (int i = 0; i < M; i++){
		Mt = lookback_option();
		if (Mt >K){
			if (St > K)Mt = Mt / St*K;
			
			double f = +(2 / (sigma*Mt*sqrt(T - t)))*normalpdf(dbs(Mt, min(St,K)) - 2.0 * a*sqrt(T - t));
			f -= 2.0 * a*pow(Mt / min(St, K), 2 * a / sigma)*normalCDF(-dbs(Mt, min(St, K))) / sigma / Mt;
			double df = +2.0*normalpdf(dbs(Mt, min(St, K)) - 2.0*a*sqrt(T - t))*(dbs(Mt, min(St, K)) - /*2.0*a*(T - t) -*/ 3.0*a*sqrt(T - t)) / (sigma*sigma*min(St, K)*Mt*(T - t));
			df += (4.0*a*a / (sigma*sigma*Mt*min(St, K)))*pow(Mt / min(St, K), 2.0*a / sigma)*normalCDF(-dbs(Mt, min(St, K)));
			res += exp(-r*(T - t))*(Mt -K)*df / f;
			/*double d1 = (log(Mt / St) - mu*(T - t)) / (sigma*(T - t));
			double d2 = (-log(Mt / St) - mu*(T - t)) / (sigma*(T - t));
			double f = (2.0 / (sigma*Mt*sqrt(T - t)))*normalpdf(d1) - (2.0*mu / (Mt*sigma*sigma))*pow(Mt / St, 2.0*mu / pow(sigma, 2))*normalCDF(d2);
			double df = (2.0 / (Mt*sigma*sigma*sqrt(T - t)*St))*normalpdf(d1)*(d1 - mu*sqrt(T - t) / sigma) + (4.0*mu*mu / (Mt*pow(sigma, 4)*St))*pow(Mt / St, 2.0*mu / pow(sigma, 2))*normalCDF(d2);
			res += exp(-r*(T - t))*(Mt - K)*df / f;*/
		}
	}
	return res / M;
}
double Lookback_option::gammaLRPW()const{
	double res = 0.0, Mt = 0.0;
	double a = (r - pow(sigma, 2) / 2.0) / sigma;
	for (int i = 0; i < M; i++){
		Mt = lookback_option();
		if (Mt > K){
			double d2 = dbs(Mt, St) - 2.0*a*sqrt(T - t), xxx = pow(Mt / St, 2.0*a / sigma);
			double f = +(2 / (sigma*Mt*sqrt(T - t)))*normalpdf(dbs(Mt, St) - 2.0 * a*sqrt(T - t));
			f -= 2.0 * a*pow(Mt / St, 2 * a / sigma)*normalCDF(-dbs(Mt, St)) / sigma / Mt;
			double df = +2.0*normalpdf(dbs(Mt, St) - 2.0*a*sqrt(T - t))*(dbs(Mt, St) - /*2.0*a*(T - t) -*/ 3.0*a*sqrt(T - t)) / (sigma*sigma*St*Mt*(T - t));
			df += (4.0*a*a / (sigma*sigma*Mt*St))*pow(Mt / St, 2.0*a / sigma)*normalCDF(-dbs(Mt, St));
			double ddf = (-2.0 / (St*St*pow(sigma, 2)*Mt*(T - t)))*normalpdf(d2)*(d2 - a*sqrt(T - t));
			ddf += (2.0 / (St*St*pow(sigma, 3)*Mt*pow(T - t, 1.5)))*normalpdf(d2)*(d2*(d2 - a*sqrt(T - t)) - 1);
			ddf -= (1 + 2.0*a / sigma - 1 / (sigma*sqrt(T - t)))*(4.0*pow(a, 2) / (pow(sigma*St, 2)*Mt))*pow(Mt / St, 2.0*a / sigma)*normalCDF(-dbs(Mt, St));
			if (Mt > K){
				res += exp(-r*(T - t))*Mt / St*df / f;
				res += exp(-r*(T - t))*(Mt - K)*(ddf / f - df*df / f / f);
			}
		}
	}
	return res / M;
}
double Lookback_option::gammaPWLR()const{
	double res = 0.0, Mt = 0.0;
	double a = (r - pow(sigma, 2) / 2.0) / sigma;
	if (St > K)return 0;
	for (int i = 0; i < M; i++){
		Mt = lookback_option();
		if (Mt > K){
			double f = +(2 / (sigma*Mt*sqrt(T - t)))*normalpdf(dbs(Mt, St) - 2.0 * a*sqrt(T - t));
			f -= 2.0 * a*pow(Mt / St, 2 * a / sigma)*normalCDF(-dbs(Mt, St)) / sigma / Mt;
			double df = +2.0*normalpdf(dbs(Mt, St) - 2.0*a*sqrt(T - t))*(dbs(Mt, St) - 2.0*a*(T - t) - a*sqrt(T - t)) / (sigma*sigma*St*Mt*(T - t));
			df += (4.0*a*a / (sigma*sigma*Mt*St))*pow(Mt / St, 2.0*a / sigma)*normalCDF(-dbs(Mt, St));

			res += exp(-r*(T - t))*(Mt/ St)*(df/f-1/St);
		}
	}
	return res / M;
}
double Lookback_option::gammaLR()const{
	double res = 0.0, Mt;
	double a = (r - pow(sigma, 2) / 2.0) / sigma;

	if (St <= K){


		for (int i = 0; i < M; i++){
			Mt = lookback_option();
			if (Mt >K){
				double f = +(2 / (sigma*Mt*sqrt(T - t)))*normalpdf(dbs(Mt, St) - 2.0 * a*sqrt(T - t));
				f -= 2.0 * a*pow(Mt / St, 2 * a / sigma)*normalCDF(-dbs(Mt, St)) / sigma / Mt;
				double d2 = dbs(Mt, St) - 2.0*a*sqrt(T - t), xxx = pow(Mt / St, 2.0*a / sigma);
				double df = (-2.0 / (St*St*pow(sigma, 2)*Mt*(T - t)))*normalpdf(d2)*(d2 - a*sqrt(T - t));
				df += (2.0 / (St*St*pow(sigma, 3)*Mt*pow(T -  t, 1.5)))*normalpdf(d2)*(d2*(d2 - a*sqrt(T - t)) - 1);
				df -= (1 + 2.0*a / sigma - 1 / (sigma*sqrt(T - t)))*(4.0*pow(a, 2) / (pow(sigma*St, 2)*Mt))*pow(Mt / St, 2.0*a / sigma)*normalCDF(-dbs(Mt, St));
				res += exp(-r*(T - t))*(Mt - K)*df / f;
			}
		}
	}
	else return 0.0;
	return res / M;
}
double Lookback_option::vegaLR()const{
	double a = (r - pow(sigma, 2) / 2.0) ;
	double res = 0.0;
	for (int i = 0; i < M; i++){
		double Mt = lookback_option();
		if (Mt > K){
			double d1 = dbs(Mt, St) - 2.0 * a*sqrt(T - t) / sigma;
			double d2 = -dbs(Mt, St);
			double dd1 = ((T - t)*(sigma*sigma + 2.0*r) - 2.0*log(Mt / St)) / (2.0*sqrt(T - t)*sigma*sigma);
			double dd2 = ((T - t)*(sigma*sigma + 2.0*r) + 2.0*log(Mt / St)) / (2.0*sqrt(T - t)*sigma*sigma);
			double f = (2 / (sigma*Mt*sqrt(T - t)))*normalpdf(d1);
			f -= 2.0 * a *pow(Mt / St, 2 * a / pow(sigma, 2))*normalCDF(d2) / (pow(sigma, 2)*Mt);
			double df= -2.0*normalpdf(d1)*((1 / (sigma*sigma*sqrt(T - t)*Mt)) + d1*dd1 / (sigma*sqrt(T - t)*Mt));
			df += 4 * r / (Mt*pow(sigma, 3))*normalCDF(d2)*pow(Mt / St, 2.0*a / pow(sigma, 2))*(1 + (2.0*a / pow(sigma, 2))*log(Mt / St));
			df -= 2.0*a / (Mt*sigma*sigma)*pow(Mt / St, 2.0*a / pow(sigma, 2))*dd2*normalpdf(d2);
			res += exp(-r*(T - t))*(Mt - K)*df / f;
			
		}
		
	}
	return res / M;
	//return 0.0;
}
double Lookback_option::vegaPW()const{
	double ST, Mt,U,res=0.0,Z;
	for (int i = 1; i < M; i++){
		
		U = rann();
		Z = rann();
		ST = St*exp((r - 0.5*sigma*sigma)*(T - t) + sigma*sqrt(T - t)*Z);
		Mt = St*exp((log(ST / St) + sqrt(pow(log(ST / St), 2) - 2.0*sigma*sigma*(T - t)*log(U))) / 2.0);
		if (Mt > K){
			double A = 2.0*((r - 0.5*pow(sigma, 2))*(T - t) + sigma*sqrt(T - t)*Z)*(-sigma*(T - t) + sqrt(T - t)*Z) - 4.0*sigma*(T - t)*log(U);
			res += Mt*0.5*(-sigma*(T - t) + sqrt(T - t)*Z + 0.5*pow(pow(log(ST / St), 2) - 2.0*sigma*sigma*(T - t)*log(U), -0.5)*A);
		}
	}
	return res / M;
	
}
double Black_Scholes::choose(int j){
	switch (j/4)
	{
	case 0:
		switch (j % 4)
		{
		case 0:return CF();
		case 1:return OP();
		default:return 0.0;
			break;
		}
	case 1:
		switch (j % 4)
		{
		case 0:return deltaCF();
		case 1:return deltaLR();
		case 2:return deltaPW();
		default:return 0.0;
			break;
		}
	case 2:
		switch (j % 4)
		{
		case 0:return gammaCF();
		case 1:return gammaLR();
		case 2:return gammaPWLR();
		case 3:return gammaLRPW();
		default:return 0.0;
			break;
		}
	case 3:
		switch (j % 4)
		{
		case 0:return vegaCF();
		case 1:return vegaLR();
		case 2:return vegaPW();
		default:return 0.0;
			break;
		}
	default:
		break;
	}
}
void Errortest(Black_Scholes &a,int j, pair<double, double>&z){
	double ermean = 0, erms = 0;
	int k = j / 4;
	//static int m = 0;
	double tem,y=a.choose(4*k);
	//ofstream fout;
	//m++;
	//fout.open("test"+to_string(m)+".csv"/*, ios::app*/);
	
	for (int i = 0; i < 50; i++){
		tem = a.choose(j);
		//fout << tem << "," << y << endl;
		ermean += tem-y;
		erms += pow(tem-y,2);
	}
	//fout.close();
	z= make_pair(ermean, erms);
}
void Black_Scholes::ETdelta(){
	pair<double, double>res1, res2;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	clock_t start = clock();
	thread t1(Errortest, *this,5, ref(res1));
	t1.join();
	thread t2(Errortest, *this,5, ref(res2));
	t2.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	clock_t end = clock();
	cout << "Delta by Likelihood method: " << deltaLR() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	start = clock();
	thread t3(Errortest, *this, 6, ref(res1));
	t3.join();
	thread t4(Errortest, *this, 6, ref(res2));
	t4.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	end = clock();
	cout << "Delta by Pathwise method: " << deltaPW() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
}

void Barrier_option::ETdelta(){
	pair<double, double>res1, res2;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	clock_t start = clock();
	thread t1(Errortest, *this, 5, ref(res1));
	t1.join();
	thread t2(Errortest, *this, 5, ref(res2));
	t2.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	clock_t end = clock();
	cout << "Delta by Likelihood method: " << deltaLR() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
}
void Black_Scholes::ET(){
	pair<double, double>res1, res2;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	clock_t start = clock();
	thread t1(Errortest, *this, 1, ref(res1));
	t1.join();
	thread t2(Errortest, *this, 1, ref(res2));
	t2.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	clock_t end = clock();
	cout << "Simulation price: " << OP() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
}
void Barrier_option::ET(){
	pair<double, double>res1, res2;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	clock_t start = clock();
	thread t1(Errortest, *this, 1, ref(res1));
	t1.join();
	thread t2(Errortest, *this, 1, ref(res2));
	t2.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	clock_t end = clock();
	cout << "Simulation price: " << OP() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
}

void Black_Scholes::ETgamma(){
	pair<double, double>res1, res2;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	clock_t start = clock();
	thread t1(Errortest, *this, 9, ref(res1));
	t1.join();
	thread t2(Errortest, *this, 9, ref(res2));
	t2.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	clock_t end = clock();
	cout << "Gamma by Likelihood method: " << gammaLR() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	start = clock();
	thread t3(Errortest, *this, 10, ref(res1));
	t3.join();
	thread t4(Errortest, *this, 10, ref(res2));
	t4.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	end = clock();
	cout << "Gamma by PW-LR method: " << gammaPWLR() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	start = clock();
	thread t5(Errortest, *this, 11, ref(res1));
	t5.join();
	thread t6(Errortest, *this, 11, ref(res2));
	t6.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	end = clock();
	cout << "Gamma by LR-PW method: " << gammaLRPW() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
}

void Barrier_option::ETgamma(){
	pair<double, double>res1, res2;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	clock_t start = clock();
	thread t1(Errortest, *this, 9, ref(res1));
	t1.join();
	thread t2(Errortest, *this, 9, ref(res2));
	t2.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	clock_t end = clock();
	cout << "Gamma by Likelihood method: " << gammaLR() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
}
void Black_Scholes::ETvega(){
	pair<double, double>res1, res2;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	clock_t start = clock();
	thread t1(Errortest, *this, 13, ref(res1));
	t1.join();
	thread t2(Errortest, *this, 13, ref(res2));
	t2.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	clock_t end = clock();
	cout << "vega by Likelihood method: " << vegaLR() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	start = clock();
	thread t3(Errortest, *this, 14, ref(res1));
	t3.join();
	thread t4(Errortest, *this, 14, ref(res2));
	t4.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	end = clock();
	cout << "vega by Pathwise method: " << vegaPW() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
}
void Barrier_option::ETvega(){
	pair<double, double>res1, res2;
	res1 = make_pair(0, 0), res2 = make_pair(0, 0);
	clock_t start = clock();
	thread t1(Errortest, *this, 13, ref(res1));
	t1.join();
	thread t2(Errortest, *this, 13, ref(res2));
	t2.join();
	res1.first += res2.first;
	res1.second += res2.second;
	res1.first /= 100;
	res1.second = res1.second / 100 - res1.first*res1.first;
	clock_t end = clock();
	cout << "vega by Likelihood method: " << vegaLR() << endl;
	cout << "Error mean: " << res1.first << endl << "variance: " << res1.second << endl;
	cout << "Time: " << (double)(end - start) / CLOCKS_PER_SEC << endl << endl;
}


double Double_Barrier::Barrier_option()const{
	double ST = St, mt = St, Mt = St;
	double b = (T - t) / N;
	for (int i = 0; i < N; i++){
		ST = STT(ST, t + (i + 1)*b, t + i*b);
		if (mt > ST)mt = ST;
		if (Mt < ST)Mt = ST;
		if (mt < L)return 0;
		if (Mt > U)return 0;
	}
	if (ST > K)return ST;
	else return 0;
}
double Double_Barrier::OP()const{
	double res = 0.0;
	double tem;
	for (int i = 0; i < M; i++){
		tem = Barrier_option();
		if (tem != 0)res += exp(-r*(T - t))*(tem - K);
	}
	return res / M;
}
void Double_Barrier::input(string x){
	string tem = x, tem2;
	double num;
	int i = 0;
	while (i < tem.size()){

		tem2 = "";
		for (; i < tem.size() && tem[i] != '='; i++){
			if (tem[i] != 32){
				tem2 += tem[i];
			}
		}
		for (; i < tem.size() && (tem[i]<48 || tem[i]>57); i++);
		if (i<x.size()){
			num = stod(tem.substr(i));
			for (; i < tem.size() && tem[i] != ','; i++);
			if (tem2 == "St"){ St = num; }
			else if (tem2 == "K"){ K = num; }
			else if (tem2 == "r"){ r = num; }
			else if (tem2 == "sigma"){ sigma = num; }
			else if (tem2 == "T"){ T = num; }
			else if (tem2 == "t"){ t = num; }
			else if (tem2 == "N"){ N = num; }
			else if (tem2 == "U"){ U = num; }
			else if (tem2 == "L"){ L = num; }
			else if (tem2 == "M"){ M = num; }
			i++;
		}

	}
}