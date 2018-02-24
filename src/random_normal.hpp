/*
 Random normal class.
 This class encapsulates all the services for random normal generation.
 */
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#define TIME_USED(f,s) ((double)((s)-(f))/(double)CLOCKS_PER_SEC)
using namespace std;

class random_normal{
    // VARIABLES
    double mean, variance;
    vector<long double> normal_random_numbers;
    long double sum = 0, sum_of_squares = 0;
    // UTILITY FUNCTIONS
    
    // MAIN FUNCTIONS
    
    // Generates m numbers from distribution N(mean,variance) using Marsaglia polar method
    void generate_by_marsaglia (long int m);
public:
    // Constructor
    random_normal(double random_variable_mean, double random_variable_variance);
    // ACCESS FUNCTIONS
    
    // Generates m numbers from distribution N(mean,variance), shows sample statistics and execution time
    void generate (long int m);
    
    // Returns the value of i-th random number
    long double operator[] (const long int i) const;
    
    // SERVICE FUNCTIONS
};







/*
// Computes the cdf of the normal distribution
double normal_cdf(double x);
// Computes the pdf of the normal distribution
double normal_pdf(double x);
*/