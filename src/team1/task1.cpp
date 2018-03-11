#include<iostream>
#include<fstream>
#include<thread>
#include"Black_Scholes1.h"

#include<string>
//#include <iomanip>
using namespace std;
//default_random_engine generator;
//normal_distribution<double> distribution(0.0, 1.0);
#define Pi 3.141592653


//double stn(){ return distribution(generator); }
//int baic_task(){
//	Black_Scholes a(100,100,0.05,0.4,1,0);
//	/*a.changeM(100000);
//	cout << "C by closed form = " << a.black_scholes_closed_form() << endl;
//	cout << "C(0,T) = "<<a.optional_price() << endl;
//	Errortest(a, &Black_Scholes::optional_price, a.black_scholes_closed_form(), cout);
//	cout << endl << "Delta(by CF) = " << a.deltaCF() << endl;
//	cout << "Delta(by LR) = " << a.deltaLR() << endl;
//	Errortest(a, &Black_Scholes::deltaLR, a.deltaCF(), cout);
//	cout << endl << "Delta(by PW) = " << a.deltaPW() << endl;
//	Errortest(a, &Black_Scholes::deltaPW, a.deltaCF(), cout);
//	cout << endl << "Gamma(by CF) = " << a.gammaCF() << endl;
//	cout << "Gamma(by LR-PW) = " << a.gammaLRPW() << endl;
//	Errortest(a, &Black_Scholes::gammaLRPW, a.gammaCF(), cout) ;
//	cout << endl << "Gamma(by PW-LR) = " << a.gammaPWLR() << endl;
//	Errortest(a, &Black_Scholes::gammaPWLR, a.gammaCF(), cout) ;
//	cout << endl << "Gamma(by LR-LR) = " << a.gammaLRLR() << endl;
//	Errortest(a, &Black_Scholes::gammaLRLR, a.gammaCF(), cout) ;
//	cout << endl << "vega(by CF) = " << a.vegaCF() << endl;
//	cout << "vega(by LR) = " << a.vegaLR() << endl;
//	Errortest(a, &Black_Scholes::vegaLR, a.vegaCF(), cout) ;
//	cout << endl << "vega(by PW) = " << a.vegaPW() << endl;
//	Errortest(a, &Black_Scholes::vegaPW, a.vegaCF(), cout) ;
//	cout << endl;*/
//	a.changeN(10000);
//	a.changeM(10000);
//	cout << "Barrier option price(closed form) = " << a.BOCF() << endl;
//	cout << "Barrier option price(simulation) = " << a.EBO() << endl;
//	Errortest(a, &Black_Scholes::EBO, 80, a.BOCF(), cout);
//	cout << endl;
//	cout << "Barrier option delta(CF) = " << a.BOdeltaCF(80) << endl;
//	cout << "Barrier option delta(LR) = " << a.BOdeltaLR(80) << endl;
//	Errortest(a, &Black_Scholes::BOdeltaLR, 80, a.BOdeltaCF(80), cout);
//	cout << "Double Barrier option price = " << a.EBO(80, 120) << endl;
//	cout << "Lookback option price = " << a.ELB() << endl;
//	return 0;
//}
//int s0_delta(){
//	ofstream fout;
//	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
//	fout.open("s0_vs_delta.csv"/*, ios::app*/);
//	fout << "S0,deltaCF,deltaPW,deltaLR" << endl;
//	for (int i = 0; i < 101; i++){
//		a.changeSt(50+2*i);
//		fout << 50 + 2 * i << "," << a.deltaCF() << "," << a.deltaPW() << "," << a.deltaLR() << endl;
//	}
//	fout.close();
//	return 0;
//}
//int s0_gamma(){
//	ofstream fout;
//	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
//	fout.open("s0_vs_gamma.csv"/*, ios::app*/);
//	fout << "S0,gammaCF,gammaLRPW,gammaPWLR,gammaLRLR" << endl;
//	for (int i = 0; i < 101; i++){
//		a.changeSt(50 + 2 * i);
//		fout << 50 + 2 * i << "," << a.gammaCF() << "," << a.gammaLRPW() << "," << a.gammaPWLR() << "," << a.gammaLRLR() << endl;
//	}
//	fout.close();
//	return 0;
//}
//int sigma_vega(){
//	ofstream fout;
//	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
//	fout.open("sigma_vs_vega.csv"/*, ios::app*/);
//	fout << "S0,vegaCF,vegaPW,vegaLR" << endl;
//	for (int i = 0; i < 101; i++){
//		a.changesigma(0.3+0.002*i);
//		fout << 0.3 + 0.002*i << "," << a.vegaCF() << "," << a.vegaPW() << "," << a.vegaLR() << endl;
//	}
//	fout.close();
//	return 0;
//}
//int Cerror_M(){
//	ofstream fout;
//	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
//	fout.open("Cerror_M.csv"/*, ios::app*/);
//	fout << "M,errormean,errorvariance" << endl;
//	for (int i = 0; i < 101; i++){
//		a.changeM(60000+1000*i);
//		fout << 60000 + 1000 * i << ",";
//		Errortest(a, &Black_Scholes::optional_price, a.black_scholes_closed_form(), fout);
//		fout << endl;
//	}
//	fout.close();
//	return 0;
//}
//int test(double(*fun)(), ofstream& os){
//	const int z = 30;
//	int num[z] = {0},t=100000;
//	double mean = 0.0, sd = 0.0;
//	double x,a=3.5;
//	double b = a *2/ z;
//	for (int j = 0; j < t; j++){
//		x = fun();
//		/*if (x > -10 && x < 10){}
//		else{ os << x << endl; }*/
//		mean += x;
//		sd += x*x;
//		for (int i = 0; i < z; i++){
//			if ((x > -a + i*b)&&(x < -a + (i + 1)*b)){
//				num[i]++;
//			}
//		}
//	}
//	for (int i = 0; i < z; i++){
//		for (int j = 0; t * j < num[i]*200; j++){
//			cout << "*";
//		}
//		cout << endl;
//	}
//	cout <<"mean: = " << mean / t << endl;
//	cout << "sd = " << sd / t << endl;
//	return 0;
//}
//int test_1(double(*fun)(), ofstream& os){
//	const int z = 30;
//	int num[z] = { 0 }, t = 100000;
//	double mean = 0.0, sd = 0.0;
//	double x, a = 3.5;
//	double b = a * 2 / z;
//	for (int j = 0; j < t; j++){
//		x = fun();
//		if (x > -10 && x < 10){}
//		else{ os << x << endl; }
//		mean += x;
//		sd += x*x;
//		for (int i = 0; i < z; i++){
//			if ((x > -a + i*b) && (x < -a + (i + 1)*b)){
//				num[i]++;
//			}
//		}
//	}
//
//	cout << "mean: = " << mean / t << "\t";
//	cout << "sd = " << sqrt(sd / t)-1 << endl;
//	return 0;
//}
//int generator_test(){
//	ofstream fou("result.txt");
//	if (!fou.is_open()){
//		cout << "Error!\n";
//	}
//	cout << "Marsaglia polar method & output two variables:\n";
//	test_1(Black_Scholes::sndb_1, fou);
//	test_1(Black_Scholes::sndb_1, fou);
//	cout << "Marsaglia polar method & only output first variable:\n";
//	test_1(Black_Scholes::sndb, fou);
//	test_1(Black_Scholes::sndb, fou);
//	cout << "BoxMuller method & output two variable:\n";
//	test_1(Black_Scholes::snd_1, fou);
//	test_1(Black_Scholes::snd_1, fou);
//	cout << "BoxMuller method & only output first variable:\n";
//	test_1(Black_Scholes::snd, fou);
//	test_1(Black_Scholes::snd, fou);
//	cout << "from random pakeage: \n";
//	test_1(stn, fou);
//	test_1(stn, fou);
//	fou.close();
//	return 0;
//}
//double tt1(pair<double, double>&x){
//	double B = 90, r = 0.05, sigma = 0.4, St = 100,T=1,t=0;
//	double mu = r - sigma*sigma / 2;
//	double f = exp(-pow(x.first - mu*(T - t), 2) / (2 * sigma*sigma*(T - t))) / (sigma*sqrt(2 * Pi*(T - t))) - pow((B / St), 2 * mu / (sigma*sigma))*exp(-pow(x.first - 2 * B - mu*(T - t), 2) / (2 * sigma*sigma*(T - t))) / (sigma*sqrt(2 * Pi*(T - t)));
//	return f;
//}
//int test(){
//	for (int i = 0; i < 1000; i++){
//		cout << tt1(make_pair(95 + 0.01*i, 100.0))<<endl;
//	}
//	return 0;
//}
//int BOpricetest(){
//	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
//	ofstream fout;
//	a.changeN(10000);
//	a.changeM(5000);
//	fout.open("Barrier_vs_BOprice.csv"/*, ios::app*/);
//	fout << "Barrier,BO,BOCF\n";
//	for (int i = 101; i < 200; i++){
//		fout << i << "," << a.EBO(i) << "," << a.BOCF(i) << endl;
//	}
//	fout.close();
//	return 0;
//}
//int BOdeltatest(){
//	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
//	ofstream fout;
//	a.changeN(10000);
//	a.changeM(500);
//	fout.open("Barrier_vs_BO_BOdelta.csv"/*, ios::app*/);
//	fout << "Barrier,BO,BOCF,BOdeltaLR,BOdeltaCF\n";
//	for (int i = 101; i < 200; i++){
//		fout << i << "," << a.EBO(i) << "," << a.BOCF(i) << "," << a.BOdeltaLR(i) << "," << a.BOdeltaCF(i) << endl;
//	}
//	fout.close();
//	return 0;
//}
////int BOgammatest(){
////	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
////	cout << a.BOgammaLR(80) << endl;
////	cout << a.BOgammaCF(80) << endl;
////	return 0;
////}
//int BOgammatest(){
//	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
//	ofstream fout;
//	a.changeN(10000);
//	a.changeM(500);
//	fout.open("Barrier_vs_BOgamma.csv"/*, ios::app*/);
//	fout << "Barrier,BOgammaCF,BOgammaLR\n";
//	for (int i = 101; i < 200; i++){
//		fout << i << "," << a.BOgammaCF(i)<< "," << a.BOgammaLR(i) << endl;
//	}
//	fout.close();
//	return 0;
//}
//int BOgammaCFtest(){
//	Black_Scholes a(100, 100, 0.05, 0.4, 1, 0);
//	ofstream fout;
//	a.changeN(10000);
//	a.changeM(5000);
//	fout.open("Barrier_vs_BOgammaCF.csv"/*, ios::app*/);
//	fout << "Barrier,BOgammaCF\n";
//	for (int i = 101; i < 200; i++){
//		fout << i << "," << a.BOgammaCF(i) <<  endl;
//	}
//	fout.close();
//	return 0;
//}
int BOvegatest(){
	Barrier_option a(100, 100, 0.05, 0.4, 1, 0);
	ofstream fout;
	a.input("N=10000,M=500,B=80");
	fout.open("test.csv"/*, ios::app*/);
	fout << ",BOvegaCF,BOvegaLR\n";
	for (int i = 101; i < 200; i++){
		
		fout << i << "," << a.vegaCF() <<","<<a.choose(13)<< endl;
	}
	fout.close();
	return 0;
}


int LBtest(){
	Lookback_option a(100, 140, 0.05, 0.4, 1, 0);
	ofstream fout;
	fout.open("Lookback_price.csv"/*, ios::app*/);
	fout << "S0,LBCF,LB_simulation\n";
	for (int i = 80; i < 180; i++){
		a.input("St="+to_string(i));
		fout << i << "," << a.CF() << "," << a.OP() << endl;
	}
	fout.close();
	return 0;
}
int LBdeltatest(){
	Lookback_option a(100, 140, 0.05, 0.4, 1, 0);
	ofstream fout;
	fout.open("Lookback_delta.csv"/*, ios::app*/);
	fout << "S0,deltaCF,deltaPW,deltaLR\n";
	for (int i = 80; i < 180; i++){
		a.input("St=" + to_string(i));
		fout << i << "," << a.deltaCF() << "," << a.deltaPW()<<","<<a.deltaLR() << endl;
	}
	fout.close();
	return 0;
}
int LBgammatest(){
	Lookback_option a(100, 140, 0.05, 0.4, 1, 0);
	ofstream fout;
	fout.open("Lookback_gamma.csv"/*, ios::app*/);
	fout << "S0,gammaCF,gammaPWLR,gammaLRPW,gammaLR\n";
	for (int i = 80; i < 180; i++){
		a.input("St=" + to_string(i));
		fout << i << "," << a.gammaCF() << "," << a.gammaPWLR() << "," << a.gammaLRPW() << "," << a.gammaLR() << endl;
	}
	fout.close();
	return 0;
}
int input(Black_Scholes* a){
	string tem;
	pair<double, double>res1, res2;
	getline(cin, tem);
	if (tem == "Europe_Call_option"){
		delete a;
		a = new Black_Scholes;
		a->show();
		return input(a);
	}
	else if (tem == "Barrier_option"){
		delete a;
		a = new Barrier_option;
		a->show();
		return input(a);
	}
	else if (tem == "Lookback_option"){
		delete a;
		a = new Lookback_option;
		a->show();
		return input(a);
	}
	else if (tem == "Double_Barrier"){
		delete a;
		a = new Double_Barrier;
		a->show();
		return input(a);
	}
	else if (tem == "delta"){
		a->delta();
	}
	else if (tem == "gamma"){
		a->gamma();
	}
	else if (tem == "vega"){
		a->vega();
	}
	else if (tem == "price"){
		a->price();
	}
	else if (tem == "show"){
		a->show();
	}
	else if (tem == "quit"){
		return 1;
	}
	else if (tem == "Error_test_delta"){
		a->ETdelta();
	}
	else if (tem == "Error_test_price"){
		a->ET();
	}
	else if (tem == "Error_test_gamma"){
		a->ETgamma();
	}
	else if (tem == "Error_test_vega"){
		a->ETvega();
	}
	else if (tem == "linear_congruential"){
		a->changeMethod(0);
	}
	else if (tem == "mt19937"){
		a->changeMethod(1);
	}
	else if (tem == "help"){
		cout << "Europe_Call_option\n" << "Barrier_option\nshow\n" << "quit\n\n";
	}
	else{
		a->input(tem);
	}
	return input(a);
}
int write_csv(string pathToFile){
	ofstream myfile;
	myfile.open(pathToFile, ios::app);
	myfile << "iterations,delta th,delta th time,delta lr,delta lr time,gamma th,gamma th time,gamma lrlr,gamma lrlr time,vega th,vega th time,vega lr,vega lr time\n";

	Barrier_option s1(100, 100, 0.05, 0.4, 1, 0 );

	for (int i = 1; i<100000; i *= 10){
		for (int j = 0; j<1000; j += 1){
			s1.input("M=" + to_string(i));
			clock_t begin = clock();
			double delta_th = s1.deltaCF();
			clock_t end = clock();
			double delta_th_time = double(end - begin) / CLOCKS_PER_SEC;
			begin = clock();
			double delta_lr = s1.deltaLR();
			end = clock();
			double delta_lr_time = double(end - begin) / CLOCKS_PER_SEC;



			myfile << i << "," << delta_th << "," << delta_th_time;
			myfile << "," << delta_lr << "," << delta_lr_time;

			begin = clock();
			double gamma_th = s1.gammaCF();
			end = clock();
			double gamma_th_time = double(end - begin) / CLOCKS_PER_SEC;
			begin = clock();
			double gamma_lrlr = s1.gammaLR();
			end = clock();
			double gamma_lrlr_time = double(end - begin) / CLOCKS_PER_SEC;


			myfile << "," << gamma_th << "," << gamma_th_time;
			myfile << "," << gamma_lrlr << "," << gamma_lrlr_time;

			begin = clock();
			double vega_th = s1.vegaCF();
			end = clock();
			double vega_th_time = double(end - begin) / CLOCKS_PER_SEC;
			begin = clock();
			double vega_lr = s1.vegaLR();
			end = clock();
			double vega_lr_time = double(end - begin) / CLOCKS_PER_SEC;
			

			myfile << "," << vega_th << "," << vega_th_time;
			myfile << "," << vega_lr << "," << vega_lr_time;
			


			myfile << "\n";

		}
	}


	return 0;
}
int main(){
	cout << (char)07;
	//generator_test();
	//s0_delta();
	//s0_gamma();
	//sigma_vega();
	//Cerror_M();
	//baic_task();
	//test();
	//BOpricetest();
	//BOdeltatest();
	//BOgammatest();
	//BOvegatest2();
	//LBtest();
	//LBdeltatest();
	//LBgammatest();
	Black_Scholes *a;
	a = new Black_Scholes;
	a->show();
	input(a);
	//write_csv("summary.csv");
	cout << (char)07;
	return 0;
}