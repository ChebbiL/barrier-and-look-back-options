#ifndef RANDOM_NORMAL_H  // header guard
#define RANDOM_NORMAL_H
#include<iostream>  // For console output
#include <vector>   // For storing random numbers
#include <cmath>    // For pow() function and other math
#include <random>   //
#include <ctime>    //
#ifndef M_PI
#define M_PI 3.1415926535898
#endif

#define TIME_USED(f,s) ((double)((s)-(f))/(double)CLOCKS_PER_SEC) // For computing implementation time

using namespace std;

// GLOBAL FUNCTIONS

// Computes the cdf of the normal distribution
double normal_cdf(double x);

// Computes the pdf of the normal distribution
double normal_pdf(double x);

// Generates a single random number from U[0,1]
double get_Urandom(void);
// Generates a single random number from N(0,1)
double get_Nrandom(void);

/*
 Random normal class.
 This class encapsulates the services for normal random variables, namely:
 1) computing the cdf/pdf;
 2) random generation
 */

class random_normal {
    // VARIABLES
    double mean, variance; // parameters of distribution, N(mean, variance)
    vector<double> normal_random_numbers;
    long double sum = 0, sum_of_squares = 0;
    double sample_mean, sample_variance;
    double time_used;
    // UTILITY FUNCTIONS

    void report (long int m, double time_used);
    // MAIN FUNCTIONS

    // Generates m numbers from distribution N(mean,variance) using Marsaglia polar method
    void generate_by_marsaglia (long int m);
public:
    // Constructor
    random_normal(double random_variable_mean = 0, double random_variable_variance = 1){
        mean = random_variable_mean; variance = random_variable_variance;
    }
    // ACCESS FUNCTIONS

    // Generates m numbers from distribution N(mean,variance), shows sample statistics and execution time
    void generate (long int m);

    // Returns the value of i-th random number
    double operator[] (const long int i) const;

    // SERVICE FUNCTIONS
};


// Computes the cdf of the normal distribution
double normal_cdf(double x);
// Computes the pdf of the normal distribution
double normal_pdf(double x);


#endif
