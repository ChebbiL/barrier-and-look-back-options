#include<iostream>
#include<fstream>
#include"Black_Scholes1.h"
#include <random>
#include <iomanip>
using namespace std;
default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.0);
#define Pi 3.141592653


double stn(){ return distribution(generator); }
int baic_task(){
	Black_Scholes a(100,100,0.05,0.4,1,0);
	//a.changeM(100000);
	//cout << "C by closed form = " << a.black_scholes_closed_form() << endl;
	//cout << "C(0,T) = "<<a.optional_price() << endl;
	//Errortest(a, &Black_Scholes::optional_price, a.black_scholes_closed_form(), cout);
	//cout << endl << "Delta(by CF) = " << a.deltaCF() << endl;
	//cout << "Delta(by LR) = " << a.deltaLR() << endl;
	//Errortest(a, &Black_Scholes::deltaLR, a.deltaCF(), cout);
	//cout << endl << "Delta(by PW) = " << a.deltaPW() << endl;
	//Errortest(a, &Black_Scholes::deltaPW, a.deltaCF(), cout);
	//cout << endl << "Gamma(by CF) = " << a.gammaCF() << endl;
	//cout << "Gamma(by LR-PW) = " << a.gammaLRPW() << endl;
	//Errortest(a, &Black_Scholes::gammaLRPW, a.gammaCF(), cout) ;
	//cout << endl << "Gamma(by PW-LR) = " << a.gammaPWLR() << endl;
	//Errortest(a, &Black_Scholes::gammaPWLR, a.gammaCF(), cout) ;
	//cout << endl << "Gamma(by LR-LR) = " << a.gammaLRLR() << endl;
	//Errortest(a, &Black_Scholes::gammaLRLR, a.gammaCF(), cout) ;
	//cout << endl << "vega(by CF) = " << a.vegaCF() << endl;
	//cout << "vega(by LR) = " << a.vegaLR() << endl;
	//Errortest(a, &Black_Scholes::vegaLR, a.vegaCF(), cout) ;
	//cout << endl << "vega(by PW) = " << a.vegaPW() << endl;
	//Errortest(a, &Black_Scholes::vegaPW, a.vegaCF(), cout) ;
	//cout << endl;
	a.changeN(1000);
	cout << "Barrier option price(simulation) = " << a.EBO(90) << endl;
	cout << "Barrier option price(closed form) = " << a.BOCF(90) << endl;
	cout << "Barrier option delta(LR) = " << a.BOdeltaLR(90) << endl;
	cout << "Barrier option delta(CF) = " << a.BOdeltaCF(90) << endl;
	cout << "Double Barrier option price = " << a.EBO(90, 120) << endl;
	cout << "Lookback option price = " << a.ELB(100) << endl;
	return 0;
}
int s0_delta(){
	ofstream fout;
	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
	fout.open("s0_vs_delta.csv"/*, ios::app*/);
	fout << "S0,deltaCF,deltaPW,deltaLR" << endl;
	for (int i = 0; i < 101; i++){
		a.changeSt(50+2*i);
		fout << 50 + 2 * i << "," << a.deltaCF() << "," << a.deltaPW() << "," << a.deltaLR() << endl;
	}
	fout.close();
	return 0;
}
int s0_gamma(){
	ofstream fout;
	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
	fout.open("s0_vs_gamma.csv"/*, ios::app*/);
	fout << "S0,gammaCF,gammaLRPW,gammaPWLR,gammaLRLR" << endl;
	for (int i = 0; i < 101; i++){
		a.changeSt(50 + 2 * i);
		fout << 50 + 2 * i << "," << a.gammaCF() << "," << a.gammaLRPW() << "," << a.gammaPWLR() << "," << a.gammaLRLR() << endl;
	}
	fout.close();
	return 0;
}
int sigma_vega(){
	ofstream fout;
	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
	fout.open("sigma_vs_vega.csv"/*, ios::app*/);
	fout << "S0,vegaCF,vegaPW,vegaLR" << endl;
	for (int i = 0; i < 101; i++){
		a.changesigma(0.3+0.002*i);
		fout << 0.3 + 0.002*i << "," << a.vegaCF() << "," << a.vegaPW() << "," << a.vegaLR() << endl;
	}
	fout.close();
	return 0;
}
int Cerror_M(){
	ofstream fout;
	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
	fout.open("Cerror_M.csv"/*, ios::app*/);
	fout << "M,errormean,errorvariance" << endl;
	for (int i = 0; i < 101; i++){
		a.changeM(60000+1000*i);
		fout << 60000 + 1000 * i << ",";
		Errortest(a, &Black_Scholes::optional_price, a.black_scholes_closed_form(), fout);
		fout << endl;
	}
	fout.close();
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
	test_1(Black_Scholes::sndb_1, fou);
	test_1(Black_Scholes::sndb_1, fou);
	cout << "Marsaglia polar method & only output first variable:\n";
	test_1(Black_Scholes::sndb, fou);
	test_1(Black_Scholes::sndb, fou);
	cout << "BoxMuller method & output two variable:\n";
	test_1(Black_Scholes::snd_1, fou);
	test_1(Black_Scholes::snd_1, fou);
	cout << "BoxMuller method & only output first variable:\n";
	test_1(Black_Scholes::snd, fou);
	test_1(Black_Scholes::snd, fou);
	cout << "from random pakeage: \n";
	test_1(stn, fou);
	test_1(stn, fou);
	fou.close();
	return 0;
}
double tt1(pair<double, double>&x){
	double B = 90, r = 0.05, sigma = 0.4, St = 100,T=1,t=0;
	double mu = r - sigma*sigma / 2;
	double f = exp(-pow(x.first - mu*(T - t), 2) / (2 * sigma*sigma*(T - t))) / (sigma*sqrt(2 * Pi*(T - t))) - pow((B / St), 2 * mu / (sigma*sigma))*exp(-pow(x.first - 2 * B - mu*(T - t), 2) / (2 * sigma*sigma*(T - t))) / (sigma*sqrt(2 * Pi*(T - t)));
	return f;
}
int test(){
	for (int i = 0; i < 1000; i++){
		cout << tt1(make_pair(95 + 0.01*i, 100.0))<<endl;
	}
	return 0;
}
int main(){
	//generator_test();
	//s0_delta();
	//s0_gamma();
	//sigma_vega();
	//Cerror_M();
	//baic_task();
	test();
	return 0;
}