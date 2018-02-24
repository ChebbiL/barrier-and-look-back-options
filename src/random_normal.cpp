#include "random_normal.hpp"


random_device rd;             // an engine used to get a SEED for Mersenne Twister
mt19937 mt(rd());             //  Mersenne Twister engine
uniform_real_distribution<long double> dist(0, 1);


long double get_random(void) {
	return dist(mt);
}

// Generates m numbers from distribution N(mean,variance) using Marsaglia polar method
void random_normal::generate_by_marsaglia (long int m) {
	normal_random_numbers.resize(m, 0);
    if (m%2==1) normal_random_numbers.resize(m+1, 0);
    long double v1, v2, w;
    for (int j = 0; j < m; j += 2) {
        w = 1.1;  // set w>1
        while (w > 1) {
            v1 = 2.0 * dist(mt) - 1;
            v2 = 2.0 * dist(mt) - 1;
            w = v1*v1 + v2*v2;
        }
        normal_random_numbers[j]      = sqrt(-2 * log(w) / w)*v1;
        normal_random_numbers[j + 1]  = sqrt(-2 * log(w) / w)*v2;
        normal_random_numbers[j]     *= sqrt(variance);
        normal_random_numbers[j + 1] *= sqrt(variance);
        normal_random_numbers[j]     += mean;
        normal_random_numbers[j+1]   += mean;
        sum            += (normal_random_numbers[j] + normal_random_numbers[j + 1]);                  //sum of normal RVs
        sum_of_squares += pow(normal_random_numbers[j], 2) + pow(normal_random_numbers[j + 1], 2);    //sum squared of normal RVs
    }
}

// Generates m numbers from distribution N(mean,variance), shows sample statistics and execution time
void random_normal::generate (long int m) {
    clock_t first = clock();      // clock starts
    generate_by_marsaglia(m);
    double sample_mean = sum / m;
    double sample_variance = (sum_of_squares / m) - pow(sample_mean, 2);
    double time_used = TIME_USED(first, clock());
    report(m, time_used);
}

// Reports sample statistics and execution time
void random_normal::report (long int m)  {
    cout << m << " random variables have been generated from normal distribution"
         << "with Variance = "<<variance<<", Mean = "<<mean<<endl;
    cout << "Here is statistics of the sample: \n";
    cout << "Mean       = " << sample_mean << endl;
    cout << "Variance   = " << sample_variance << endl;
    cout << "Time used  = " << time_used << " s." << endl<<endl;     //clock stops
}

// Overloading [] operator for easy output of random numbers
long double random_normal::operator[] (const long int i) const {
    if (normal_random_numbers.size()>i) return normal_random_numbers[i];
    else {
        cout<<"Error. Vector size is less than the value entered."<<endl;
        return 0;
    }
}

<<<<<<< HEAD
double normal_cdf(double x){
  return erfc( - x / sqrt(2)) / 2;
}
double normal_pdf(double x){
  return exp( - x * x / 2) / sqrt(2 * M_PI);
}
//temporary dummy function
double fun(){
  return .5;
}
=======
// Computes the cdf of the normal distribution
double normal_cdf(double x){
  return erfc( - x / sqrt(2)) / 2;
}

// Computes the cdf of the normal distribution
double normal_pdf(double x){
  return exp( - x * x / 2) / sqrt(2 * M_PI);
}
>>>>>>> 4a6688d0fe0e9d2fa658064a655bfb019af730e1
