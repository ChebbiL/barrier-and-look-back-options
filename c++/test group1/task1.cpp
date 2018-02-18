#include<iostream>
#include<fstream>
#include"Black_Scholes1.h"
#include <random>
using namespace std;
default_random_engine generator;
normal_distribution<double> distribution(0.0, 1.0);
#define Pi 3.141592653


double stn(){ return distribution(generator); }
int baic_task(){
	Black_Scholes a(100,100,0.05,0.4,1,0);
	a.changeM(10000);
	cout << "C by closed form = " << a.black_scholes_closed_form() << endl;
	cout << "C(0,T) = "<<a.optional_price() << endl;
	Errortest(a, &Black_Scholes::optional_price, a.black_scholes_closed_form());
	cout << "Delta(by CF) = " << a.deltaCF() << endl;
	cout << "Delta(by LR) = " << a.deltaLR() << endl;
	Errortest(a, &Black_Scholes::deltaLR, a.deltaCF());
	cout << "Delta(by PW) = " << a.deltaPW() << endl;
	Errortest(a, &Black_Scholes::deltaPW, a.deltaCF());
	cout << "Gamma(by CF) = " << a.gammaCF() << endl;
	cout << "Gamma(by LR-PW) = " << a.gammaLRPW() << endl;
	Errortest(a, &Black_Scholes::gammaLRPW, a.gammaCF());
	cout << "Gamma(by PW-LR) = " << a.gammaPWLR() << endl;
	Errortest(a, &Black_Scholes::gammaPWLR, a.gammaCF());
	cout << "Gamma(by LR-LR) = " << a.gammaLRLR() << endl;
	Errortest(a, &Black_Scholes::gammaLRLR, a.gammaCF());
	cout << "vega(by CF) = " << a.vegaCF() << endl;
	cout << "vega(by LR) = " << a.vegaLR() << endl;
	Errortest(a, &Black_Scholes::vegaLR, a.vegaCF());
	cout << "vega(by PW) = " << a.vegaPW() << endl;           
	Errortest(a, &Black_Scholes::vegaPW, a.vegaCF());
	
	cout << "Barrier option price = " << a.EBO(90) << endl;
	cout << "Double Barrier option price = " << a.EBO(90, 120) << endl;
	cout << "Lookback option price = " << a.ELB(100) << endl;
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


int main(){
	//generator_test();
	baic_task();
	
	return 0;
}