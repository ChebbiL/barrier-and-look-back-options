#include<iostream>
#include<fstream>
#include"Black_Scholes.h"
#include <random>
using namespace std;
default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.0);
#define Pi 3.141592653


double stn(){ return distribution(generator); }
int baic_task(double(*fun)()){
	cout << "C(0,T) = "<<optional_price(100, 100, 0.05, 0.4, 1, 0, 100000,fun) << endl;
	cout << "Delta(by LR) = " << deltaLR(100, 100, 0.05, 0.4, 1, 0, 100000, fun) << endl;
	cout << "Delta(by PW) = " << deltaPW(100, 100, 0.05, 0.4, 1, 0, 100000, fun) << endl;
	cout << "Gamma(by LR-PW) = " << gammaLRPW(100, 100, 0.05, 0.4, 1, 0, 100000, fun) << endl;
	cout << "Gamma(by PW-LR) = " << gammaPWLR(100, 100, 0.05, 0.4, 1, 0, 100000, fun) << endl;
	cout << "Gamma(by LR-LR) = " << gammaLRLR(100, 100, 0.05, 0.4, 1, 0, 100000, fun) << endl;
	cout << "vega(by LR) = " << vegaLR(100, 100, 0.05, 0.4, 1, 0, 100000, fun) << endl;
	cout << "vega(by PW) = " << vegaPW(100, 100, 0.05, 0.4, 1, 0, 100000, fun) << endl;
	return 0;
}


int test(double(*fun)(), ofstream& os){
	const int z = 30;
	int num[z] = {0},t=100000;
	double mean = 0.0, sd = 0.0;
	double x,a=3.5;
	double b = a *2/ z;
	for (int j = 0; j < t; j++){
		x = fun();
		/*if (x > -10 && x < 10){}
		else{ os << x << endl; }*/
		mean += x;
		sd += x*x;
		for (int i = 0; i < z; i++){
			if ((x > -a + i*b)&&(x < -a + (i + 1)*b)){
				num[i]++;
			}
		}
	}
	for (int i = 0; i < z; i++){
		for (int j = 0; t * j < num[i]*200; j++){
			cout << "*";
		}
		cout << endl;
	}
	cout <<"mean: = " << mean / t << endl;
	cout << "sd = " << sd / t << endl;
	return 0;
}
int test_1(double(*fun)(), ofstream& os){
	const int z = 30;
	int num[z] = { 0 }, t = 100000;
	double mean = 0.0, sd = 0.0;
	double x, a = 3.5;
	double b = a * 2 / z;
	for (int j = 0; j < t; j++){
		x = fun();
		if (x > -10 && x < 10){}
		else{ os << x << endl; }
		mean += x;
		sd += x*x;
		for (int i = 0; i < z; i++){
			if ((x > -a + i*b) && (x < -a + (i + 1)*b)){
				num[i]++;
			}
		}
	}

	cout << "mean: = " << mean / t << "\t";
	cout << "sd = " << sqrt(sd / t)-1 << endl;
	return 0;
}
int generator_test(){
	ofstream fou("result.txt");
	if (!fou.is_open()){
		cout << "Error!\n";
	}
	cout << "Marsaglia polar method & output two variables:\n";
	test_1(sndb_1,fou);
	test_1(sndb_1, fou);
	cout << "Marsaglia polar method & only output first variable:\n";
	test_1(sndb, fou);
	test_1(sndb, fou);
	cout << "BoxMuller method & output two variable:\n";
	test_1(snd_1, fou);
	test_1(snd_1, fou);
	cout << "BoxMuller method & only output first variable:\n";
	test_1(snd, fou);
	test_1(snd, fou);
	cout << "from random pakeage: \n";
	test_1(stn, fou);
	test_1(stn, fou);
	fou.close();
	return 0;
}
//barrier B, number of partitions N
double Barrier_option(double St, double K, double r, double sigma, double T, double t, double(*fun)(), double B, int N){
	double ST = St,mt=St;
	double b = (T - t) / N;
	for (int i = 0; i < N ; i++){
		ST = STT(St, K, r, sigma, t + (i + 1)*b, t + i*b, fun);
		if (mt > ST)mt = ST;
		if (mt < B)return 0;
	}
	if (ST > K)return exp(-r*(T - t))*(ST - K);
	else return 0;
}
//barrier B, number of partitions N
double EBO(double St, double K, double r, double sigma, double T, double t, int M, double(*fun)(), double B, int N){
	double res = 0.0;
	for (int i = 0; i < M; i++){
		res += Barrier_option(St, K, r, sigma, T, t, fun, B, N);
	}
	return res / M;
}
//B1 lower bound, B2 upper bound
double Barrier_option(double St, double K, double r, double sigma, double T, double t, double(*fun)(), double B1,double B2, int N){
	double ST = St, mt = St,Mt=St;
	double b = (T - t) / N;
	for (int i = 0; i < N; i++){
		ST = STT(St, K, r, sigma, t + (i + 1)*b, t + i*b, fun);
		if (mt > ST)mt = ST;
		if (Mt < ST)Mt = ST;
		if (mt < B1)return 0;
		if (Mt > B2)return 0;
	}
	if (ST > K)return exp(-r*(T - t))*(ST - K);
	else return 0;
}
double EBO(double St, double K, double r, double sigma, double T, double t, int M, double(*fun)(), double B1,double B2, int N){
	double res = 0.0;
	for (int i = 0; i < M; i++){
		res += Barrier_option(St, K, r, sigma, T, t, fun, B1,B2, N);
	}
	return res / M;
}
double lookback_option(double St, double K, double r, double sigma, double T, double t, double(*fun)(), int N){
	double ST = St, Mt = St;
	double b = (T - t) / N;
	for (int i = 0; i < N; i++){
		ST = STT(St, K, r, sigma, t + (i + 1)*b, t + i*b, fun);
		if (Mt < ST)Mt = ST;
	}
	if (Mt > K)return exp(-r*(T - t))*(Mt - K);
	else return 0;
}
double ELB(double St, double K, double r, double sigma, double T, double t, int M, double(*fun)(), int N){
	double res = 0.0;
	for (int i = 0; i < M; i++){
		res += lookback_option(St, K, r, sigma, T, t, fun, N);
	}
	return res / M;
}
int BO_task(double(*fun)()){
	cout << "Barrier option price = " << EBO(100, 100, 0.05, 0.4, 1, 0, 100000, fun, 90, 100) << endl;
	cout << "Double Barrier option price = " << EBO(100, 100, 0.05, 0.4, 1, 0, 100000, fun, 90,120, 100) << endl;
	cout << "Lookback option price = " << ELB(100, 100, 0.05, 0.4, 1, 0, 100000, fun, 100) << endl;
	return 0;
}

int main(){
	//generator_test();
	//baic_task(snd_1);
	BO_task(snd_1);
	return 0;
}