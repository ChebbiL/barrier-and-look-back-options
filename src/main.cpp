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


int debug_thomas(){
  european_option call(100, 100, 0.05, 0.4, 1, 0, 1000);
  cout << "Price: " << call.price() << endl;
  cout << "Delta: " << call.delta() << " (LR)" << endl;
  cout << "Delta: " << call.delta("pw") << " (PW)" << endl;
  cout << "Gamma: " << call.gamma() << " (PWLR)" << endl;
  cout << "Gamma: " << call.delta("lrpw") << " (LRPW)" << endl;
  cout << "Gamma: " << call.delta("lrlr") << " (LRLR)" << endl;

  return 0;
}

int debug_konstantin(){

  return 0;
}

int main(){
  cout << "hello" << endl;
  debug_thomas();
  return 0;
}
