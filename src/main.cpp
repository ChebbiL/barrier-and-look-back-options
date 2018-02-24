/*
This is the main project file.
*/

#include<iostream> // For console output

#include "european_option.hpp"
#include "barrier_option.hpp"

using namespace std;


/*
Random normal class.
This class encapsulates all the services for random normal generation.
*/
class random_normal{
  // VARIABLES
  double mean, variance;
  // UTILITY FUNCTIONS

  // MAIN FUNCTIONS

public:
  // Constructor
  random_normal(double random_variable_mean = 0, double random_variable_variance = 1){
    mean = random_variable_mean;
    variance = random_variable_variance;
  }
  // ACCESS FUNCTIONS

  // SERVICE FUNCTIONS

};



int main(){
  cout << "hello" << endl;
  return 0;
}
