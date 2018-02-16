

unsigned long int ran();//get a large random number 2^31-1
long double rann();//get a uniform random number
double snd();   //Marsaglia polar method & using one output only
double snd_1();//using both outputs
double sndb();  //BoxMuller method & using one output only
double sndb_1();//using both outputs
double STT(double St, double K, double r, double sigma, double T, double t, double(*fun)());
double BSCall(double St, double K, double r, double sigma, double T, double t, double(*fun)());
//initial price St, strike price K, interest rate r,volatility sigma, maturity time T, starting time t
double optional_price(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)());
//repeating times M
double deltaPW(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)());
double deltaLR(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)());
double gammaLRPW(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)());
double gammaPWLR(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)());
double gammaLRLR(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)());
double vegaLR(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)());
double vegaPW(double St, double K, double r, double sigma, double T, double t, long int M, double(*fun)());

