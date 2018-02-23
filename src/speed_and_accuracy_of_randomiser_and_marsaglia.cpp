#include <iostream>
#include <ctime>
#include <random>
#define TIME_USED(f,s) ((double)((s)-(f))/(double)CLOCKS_PER_SEC)
#define SIZE 10000000 // number of normal N(0,1) variables generated
using namespace std;

int main() {
	random_device rd;             // an engine used to get a SEED for Mersenne Twister
	mt19937 mt(rd());             //  Mersenne Twister engine
	uniform_real_distribution<long double> dist(0, 1);
	//////////////////////////////////////////////////////////////////////// Checking U[0,1] random generator speed:
	clock_t first = clock();      //clock starts
	long double ld;
	for (int i = 0; i < 10000000; ++i)   ld = dist(mt);    //Generating 10MLN numbers from U[0,1]
	cout << "Time used to generate 10MLN numbers from U[0,1] : " << TIME_USED(first, clock()) << " Sec."; // clock stops
	cout << endl << endl;
	// Results:
	//   1.11 sec to generate 10 MLN N[0,1] numbers
	//////////////////////////////////////////////////////////////////////// Checking Marsaglia polar method:
	first = clock();                   //clock starts
	vector<long double> unirand(SIZE); //vector filled with normal RVs
	long double nrand1, nrand2, nsum = 0, nsumsqr = 0, v1, v2, w;
	
	for (int j = 0; j <SIZE; j += 2) {
		w = 1.1;  // set w>1
		while (w > 1) {
			v1 = 2.0 * dist(mt)-1;
			v2 = 2.0 * dist(mt)-1;
			w = v1*v1 + v2*v2;
		}
		unirand[j] = sqrt(-2 * log(w) / w)*v1;
		unirand[j+1] = sqrt(-2 * log(w) / w)*v2;

		nsum += (unirand[j] + unirand[j+1]);                                  //sum of normal RVs
		nsumsqr += (unirand[j] * unirand[j]) + (unirand[j+1] * unirand[j+1]); //sum squared of normal RVs
	}

	long double mean, var;
	mean = (nsum / SIZE);
	var = (nsumsqr / SIZE) - (mean*mean);
	cout << "Mean: " << mean << endl;
	cout << "Variance: " << var << endl; 
	cout << "Time used to generate N(0,1) RVs using Marsaglia polar method :" << TIME_USED(first, clock()) << " Sec.";  // clock stops
	cout << endl << endl;
	// Results:
	//   3.55 sec to generate 10 MLN N[0,1] numbers
	//   Mean: -0.000124
	//   Var :  0.9998



}