
 
class Black_Scholes{
private:
	double St,K,r,sigma,T,t;
    double(*fun)(void);
	int M, N;
	double BCF, deltaCF1, gammaCF1, vegaCF1;
	static unsigned long int ran();//get a large random number 2^31-1
	static long double rann();//get a uniform random number
	double normalCDF(double x);
	double normalpdf(double x);
	void getCF(){
		double d = (log(St / K) + (r + sigma*sigma / 2)*(T - t)) / (sigma*sqrt(T - t));
		BCF = normalCDF(d)*St - normalCDF(d - sigma*sqrt(T - t))*K*exp(-r*(T - t));
		deltaCF1 = normalCDF(d);
		gammaCF1 = normalpdf(d) / (St*sigma*sqrt(T - t));
		vegaCF1 = St*normalpdf(d)*sqrt(T - t);
	}
public:
	static double snd();   //Marsaglia polar method
	static double snd_1();
	static double sndb();  //BoxMuller method
	static double sndb_1();
	void changeM(int L){ M = L; }
	double STT(double St, double T, double t);
	Black_Scholes(double a = 100.0, double b = 100.0, double c=0.05,double d=0.4,double e=1,double f=0){
		St = a;
		K = b;
		r = c;
		sigma = d;
		T = e;
		t = f;
		fun = sndb_1;
		M = 100000;
		N = 100;
		getCF();
	};
	double BSCall();
	//initial price St, strike price K, interest rate r,volatility sigma, maturity time T, starting time t
	double black_scholes_closed_form()const{ return BCF; }
	double optional_price();
	//repeating times M
	double deltaCF()const{ return deltaCF1; }
	double deltaPW();
	double deltaLR();
	double gammaCF()const{ return gammaCF1; }
	double gammaLRPW();
	double gammaPWLR();
	double gammaLRLR();
	double vegaCF()const{ return vegaCF1; }
	double vegaLR();
	double vegaPW();
	double Barrier_option(double B);
	double Barrier_option(double B1, double B2);
	double EBO(double B);
	double EBO(double B1, double B2);
	double lookback_option(double L);
	double ELB(double L);
	
	
};
void Errortest(Black_Scholes x, double(Black_Scholes::*foo)(), double y);