

class Black_Scholes{
protected:
	double St,K,r,sigma,T,t;
	
	int M,type;
	unsigned long int ran()const;//get a large random number 2^31-1
	long double rann()const;//get a uniform random number
	double normalCDF(double x)const;
	double normalpdf(double x)const;
	double BSC(double S0, double k)const;
	double BSCdelta(double S0, double k)const;
	double BSCgamma(double S0, double k)const;
	double BSCvega(double S0, double k)const;
	double dbs(double x, double y)const;
	double snd()const;   //Marsaglia polar method one output only
	double snd_1()const;//Marsaglia polar method
	double sndb()const;  //BoxMuller method one output only
	double sndb_1()const;//BoxMuller method
	double fun()const{ return snd_1(); }
public:
	
	double STT(double St, double T, double t)const;
	Black_Scholes(double a = 100.0, double b = 100.0, double c=0.05,double d=0.4,double e=1,double f=0){
		St = a;
		K = b;
		r = c;
		sigma = d;
		T = e;
		t = f;
		M = 100000;
		type = 0;
	};
	virtual ~Black_Scholes();
	virtual void input(std::string x);
	virtual void show()const{
		std::cout << std::endl;
		std::cout << "Europe Call option\nSt = " << St << std::endl << "K = " << K << std::endl << "r = " << r << std::endl << "sigma = " << sigma << std::endl;
		std::cout << "T = " << T << std::endl << "t = " << t << std::endl << "M = " << M << std::endl;
		std::cout << std::endl;
	}
	virtual void price()const{
		std::cout << "price closed form result: \t \t" << CF() << std::endl << "price simulation result: \t \t" << OP() << std::endl << std::endl;
	}
	virtual void delta()const{
		std::cout << "delta closed form result: \t \t" << deltaCF() << std::endl << "delta simulation result(by PW): \t" << deltaPW() << std::endl << "delta simulation result(by LR): \t" << deltaLR() << std::endl << std::endl;
	}
	virtual void gamma()const{
		std::cout << "gamma closed form result: \t \t" << gammaCF() << std::endl << "gamma simulation result(by LRPW): \t" << gammaLRPW() << std::endl << "gamma simulation result(by PWLR): \t" << gammaPWLR() << std::endl << "gamma simulation result(by LRLR): \t" << gammaLR() << std::endl << std::endl;
	}
	virtual void vega()const{
		std::cout << "vega closed form result: \t \t" << vegaCF() << std::endl << "vega simulation result(by PW): \t \t" << vegaPW() << std::endl << "vega simulation result(by LR): \t \t" << vegaLR() << std::endl << std::endl;
	}
	//initial price St, strike price K, interest rate r,volatility sigma, maturity time T, starting time t
	virtual double CF()const{
		return BSC(St, K);
	}
	virtual double OP()const;
	//repeating times M
	virtual double deltaCF()const{ return  BSCdelta(St, K); }
	double deltaPW()const;
	virtual double deltaLR()const;
	virtual double gammaCF()const{ return BSCgamma(St, K); }
	double gammaLRPW()const;
	double gammaPWLR()const;
	virtual double gammaLR()const;
	virtual double vegaCF()const{ return BSCvega(St, K); }
	virtual double vegaLR()const;
	double vegaPW()const;
	double choose(int j);
	virtual void ET();
	virtual void ETdelta();
	virtual void ETgamma();
	virtual void ETvega();
	
};
class Barrier_option:public Black_Scholes{
private:
	double B;
	int N;
	double Barrier_option_down()const;
	double Barrier_option_up()const;
	double BOdLR(double x)const;
	double BOvLR(double x)const;
	double BOgLR(double x)const;
	double BOCF_down()const;
	double BOCF_up()const;
public:
	Barrier_option(double a = 100.0, double b = 100.0, double c = 0.05, double d = 0.4, double e = 1, double f = 0, double g = 80) :Black_Scholes(a, b, c, d, e, f), B(g){ N = 100; type = 1; }
	~Barrier_option(){}
	void input(std::string x);
	void show()const{
		std::cout << std::endl;
		std::cout << "Barrier out Call option\nSt = " << St << std::endl << "K = " << K << std::endl << "B = " << B << std::endl << "r = " << r << std::endl << "sigma = " << sigma << std::endl;
		std::cout << "T = " << T << std::endl << "t = " << t << std::endl << "M = " << M << std::endl << "N = " << N << std::endl;
		std::cout << std::endl;
	}
	void price()const{
		std::cout << "price closed form result: \t \t" << CF() << std::endl << "price simulation result: \t \t" << OP() << std::endl << std::endl;
	}
	void delta()const{
		std::cout << "delta closed form result: \t \t" << deltaCF() << std::endl << "delta simulation result(by LR): \t" << deltaLR() << std::endl << std::endl;
	}
	void gamma()const{
		std::cout << "gamma closed form result: \t \t" << gammaCF() << std::endl << "gamma simulation result(by LR): \t" << gammaLR() << std::endl << std::endl;
	}
	void vega()const{
		std::cout << "vega closed form result: \t \t" << vegaCF() << std::endl << "vega simulation result(by LR): \t \t" << vegaLR() << std::endl << std::endl;
	}
	double OP()const;
	double CF()const;
	double deltaLR()const;
	double deltaCF()const;
	double gammaLR()const;
	double gammaCF()const;
	double vegaLR()const;
	double vegaCF()const;
	void ET();
	void ETdelta();
	void ETgamma();
	void ETvega();
};
class Lookback_option :public Black_Scholes{
	int N;
public:
	Lookback_option(){}
	~Lookback_option(){}
	double lookback_option();
	double ELB();
	double LBCF();
	//double LBdLR(double M);
	/*double LBdeltaLR();
	double LBpdf(double x);*/
};
void Errortest(Black_Scholes &a, int j, std::pair<double, double>&z);



