#include <iostream>
#include <ctime>
#include <random>
#define TIME_USED(f,s) ((double)((s)-(f))/(double)CLOCKS_PER_SEC)
#define PI 3.14159265359
#define M 1000000                        // number of MC simulations, can be changed
using namespace std;

void generate_N_RV();                    // uses Marsaglia polar method to generate N(0,1) RVs
void Black_Scholes_Theory_t_0(double&s0, double&strike, double &rate, double &vol, double &T);  // calculates theoretical values of call price and greeks at time 0
void Monte_Carlo(double&s0, double&strike, double &rate, double &vol, double &T);               // calculates call price and greeks at time 0 using Monte-Carlo simulation

vector<long double> unirand(M);          // vector to be filled with normal RVs

   //////// the following can be changed:
double s0 = 100;                         // initial asset price
double strike = 100;                     // strike price
double rate = 0.05;                      // interest rate
double vol = 0.4;                        // volatility
double T = 1;                            // maturity time




int main() {
	cout << fixed;
	Black_Scholes_Theory_t_0(s0, strike, rate, vol, T);   // calculates theoretical values of call price and greeks at time 0
	generate_N_RV();                                      // uses Marsaglia polar method to generate N(0,1) RVs
	Monte_Carlo(s0, strike, rate, vol, T);				  // calculates call price and greeks at time 0 using Monte-Carlo simulation
	
}







void generate_N_RV() {            // using Marsaglia polar method to generate N(0,1) RVs
	clock_t first = clock();      // clock starts
	random_device rd;             // an engine used to get a SEED for Mersenne Twister
	mt19937 mt(rd());             //  Mersenne Twister engine
	uniform_real_distribution<long double> dist(0, 1);

	long double nsum = 0, nsumsqr = 0, v1, v2, w;

	for (int j = 0; j <M; j += 2) {
		w = 1.1;  // set w>1
		while (w > 1) {
			v1 = 2.0 * dist(mt) - 1;
			v2 = 2.0 * dist(mt) - 1;
			w = v1*v1 + v2*v2;
		}
		unirand[j] = sqrt(-2 * log(w) / w)*v1;
		unirand[j + 1] = sqrt(-2 * log(w) / w)*v2;

		nsum    += (unirand[j] + unirand[j + 1]);                                    //sum of normal RVs
		nsumsqr += (unirand[j] * unirand[j]) + (unirand[j + 1] * unirand[j + 1]);    //sum squared of normal RVs
	}

	long double mean, var;
	mean = (nsum / M);
	var = (nsumsqr / M) - (mean*mean);
	cout << M << " normal random variables have been generated.\nHere is statistics of the sample: \n";
	cout << "Mean       = " << mean << endl;
	cout << "Variance   = " << var << endl;
	cout << "Time used  = " << TIME_USED(first, clock()) << " s." << endl;           //clock stops
	cout << endl << endl;
}



void Black_Scholes_Theory_t_0(double&s0, double&strike, double &rate, double &vol, double &T){   // calculates theoretical values of call price and greeks at time 0
	double call_t, delta_t, gamma_t, vega_t;  // theoretical values of call price and greeks at time 0
	double d1, d2, phi_d1, phi_d2;
	cout << endl;
	d1 = (log(s0/strike)+(rate+0.5*pow(vol,2))*T)/vol*sqrt(T);
	d2 = d1 - vol*sqrt(T);
	phi_d1 = erfc(-d1 / sqrt(2)) / 2;
	phi_d2 = erfc(-d2 / sqrt(2)) / 2;
	call_t = s0*phi_d1 - strike*exp(-rate*T)*phi_d2;
	delta_t = phi_d1;
	gamma_t = (exp(-0.5*d1*d1) / sqrt(2 * PI)) / (s0*vol*sqrt(T));
	vega_t = s0*sqrt(T)*(exp(-0.5*d1*d1) / sqrt(2 * PI));

	cout <<"Results using Black-Scholes closed-form formulas:\nCall Price = "<< call_t << endl;
	cout << "Delta      = " << delta_t << endl;
	cout << "Gamma      = " <<gamma_t << endl;
	cout << " Vega      = "<< vega_t << endl<<endl;
}




void Monte_Carlo(double&s0, double&strike, double &rate, double &vol, double &T) {    // calculates call price and greeks at time 0 using Monte-Carlo simulation
	clock_t first = clock();      //clock starts
	double sT, priceSum = 0, priceSumSqr = 0;
	double LRdeltaSum = 0, LRdeltaSumSqr = 0, PWdeltaSum = 0, PWdeltaSumSqr = 0, LRvegaSum = 0, LRvegaSumSqr = 0, PWvegaSum = 0, PWvegaSumSqr = 0;
	double LRgammaSum = 0, LRgammaSumSqr = 0, LR_PWgammaSum = 0, LR_PWgammaSumSqr = 0, PW_LRgammaSum = 0, PW_LRgammaSumSqr = 0;
	for (int j = 0; j < M; j++) {
		sT = s0*exp(rate - 0.5*vol*vol + vol*sqrt(T)*unirand[j]);                                                                       // change
		if (sT - strike > 0) {                    // if the option is exercised at maturity
			priceSum += ((sT - strike)*exp(-rate*T));
			priceSumSqr += ((sT - strike)*exp(-rate*T)) * ((sT - strike)*exp(-rate*T));
			LRdeltaSum += ((sT - strike)*exp(-rate*T))*unirand[j] / (s0*vol*sqrt(T));
			LRdeltaSumSqr += (((sT - strike)*exp(-rate*T))*unirand[j] / (s0*vol*sqrt(T)))*(((sT - strike)*exp(-rate*T))*unirand[j] / (s0*vol*sqrt(T)));
			PWdeltaSum += exp(-rate*T)*sT / s0;
			PWdeltaSumSqr += (exp(-rate*T)*sT / s0)*(exp(-rate*T)*sT / s0);
			LRgammaSum += ((sT - strike)*exp(-rate*T))* ((unirand[j] * unirand[j] - 1) / (s0*s0*vol*vol*T) - unirand[j] / (s0*s0*vol*sqrt(T)));
			LRgammaSumSqr += pow(((sT - strike)*exp(-rate*T))* ((unirand[j] * unirand[j] - 1) / (s0*s0*vol*vol*T) - unirand[j] / (s0*s0*vol*sqrt(T))), 2);
			LR_PWgammaSum += exp(-rate*T)*strike*unirand[j] / (s0*s0*vol*sqrt(T));
			LR_PWgammaSumSqr += pow((exp(-rate*T)*strike*unirand[j] / (s0*s0*vol*sqrt(T))), 2);
			PW_LRgammaSum += exp(-rate*T)*sT * (unirand[j] / (vol*sqrt(T)) - 1) / (s0*s0);
			PW_LRgammaSumSqr += pow((exp(-rate*T)*sT * (unirand[j] / (vol*sqrt(T)) - 1) / (s0*s0)), 2);
			LRvegaSum += ((sT - strike)*exp(-rate*T)* ((unirand[j] * unirand[j] - 1) / vol - unirand[j] * sqrt(T)));
			LRvegaSumSqr += pow(((sT - strike)*exp(-rate*T)* ((unirand[j] * unirand[j] - 1) / vol - unirand[j] * sqrt(T))), 2);
			PWvegaSum += exp(-rate*T)*sT*  (-vol*T + sqrt(T)*unirand[j]);
			PWvegaSumSqr += pow((exp(-rate*T)*sT*  (-vol*T + sqrt(T)*unirand[j])), 2);
		}
	}
	double priceMean, priceVar, LRdeltaMean, LRdeltaVar, PWdeltaMean, PWdeltaVar, LRvegaMean, LRvegaVar, PWvegaMean, PWvegaVar;
	double LRgammaMean, LRgammaVar, LR_PWgammaMean, LR_PWgammaVar, PW_LRgammaMean, PW_LRgammaVar;
	priceMean = priceSum / M;                            // average discounted payoff at time 0 (MC estimation of call price)
	priceVar = (priceSumSqr / M) - priceMean*priceMean;  // variance of call price
	LRdeltaMean = LRdeltaSum / M;
	LRdeltaVar = (LRdeltaSumSqr / M) - LRdeltaMean*LRdeltaMean;
	PWdeltaMean = PWdeltaSum / M;
	PWdeltaVar = (PWdeltaSumSqr / M) - PWdeltaMean*PWdeltaMean;
	LRgammaMean = LRgammaSum / M;
	LRgammaVar = (LRgammaSumSqr / M) - LRgammaMean*LRgammaMean;
	LR_PWgammaMean = LR_PWgammaSum / M;
	LR_PWgammaVar = (LR_PWgammaSumSqr / M) - LR_PWgammaMean * LR_PWgammaMean;
	PW_LRgammaMean = PW_LRgammaSum / M;
	PW_LRgammaVar = (PW_LRgammaSumSqr / M) - PW_LRgammaMean * PW_LRgammaMean;
	LRvegaMean = LRvegaSum / M;
	LRvegaVar = (LRvegaSumSqr / M) - LRvegaMean*LRvegaMean;
	PWvegaMean = PWvegaSum / M;
	PWvegaVar = (PWvegaSumSqr / M) - PWvegaMean*PWvegaMean;

	cout << "Call Price mean (MC call price estimate)    = " << priceMean << endl;
	cout << "Call Price variance                         = " << priceVar << endl;
	cout << "LR delta mean     = " << LRdeltaMean << endl;
	cout << "PW delta mean     = " << PWdeltaMean << endl;
	cout << "LR delta variance    = " << LRdeltaVar << endl;
	cout << "PW delta variance    = " << PWdeltaVar << endl;
	cout << "LR gamma mean     = " << LRgammaMean << endl;
	cout << "LR_PW gamma mean  = " << LR_PWgammaMean << endl;
	cout << "PW_LR gamma mean  = " << PW_LRgammaMean << endl;
	cout << "LR gamma variance    = " << LRgammaVar << endl;
	cout << "LR_PW gamma variance = " << LR_PWgammaVar << endl;
	cout << "PW_LR gamma variance = " << PW_LRgammaVar << endl;
	cout << "LR vega mean      = " << LRvegaMean << endl;
	cout << "PW vega mean      = " << PWvegaMean << endl;
	cout << "LR vega variance     = " << LRvegaVar << endl;
	cout << "PW vega variance     = " << PWvegaVar << endl;
	cout << "Time used         = " << TIME_USED(first, clock()) << " s.";  //clock stops
	cout << endl << endl;
}