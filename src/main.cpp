/*
This is the main project file.
*/

#include<iostream> // For console output


#include "european_option.hpp"
#include "barrier_option.hpp"


#include <cmath> // For math calculations
#include <string> // For user instructions in greeks
#include <random> // For random numbers generation

using namespace std;





int debug_thomas(int x){
  cout << "NUMBER OF SIMS: " << x << endl;
  european_option call(100, 100, 0.05, 0.4, 1, 0, x);
  cout<<"Theoretic price: "<<call.payoff_theoretic()<<endl;
  cout<<"Theoretic delta: "<<call.delta("th")<<endl;
  cout<<"Theoretic gamma: "<<call.gamma("th")<<endl;
  cout<<" Theoretic vega: "<<call.vega("th")<<endl;

  cout << "Simulation:"<<endl;
  cout << "Price: " << call.price() << endl;
  cout << "Delta: " << call.delta() << " (LR)" << endl;
  cout << "Delta: " << call.delta("pw") << " (PW)" << endl;
  cout << "Gamma: " << call.gamma() << " (PWLR)" << endl;
  cout << "Gamma: " << call.gamma("lrpw") << " (LRPW)" << endl;
  cout << "Gamma: " << call.gamma("lrlr") << " (LRLR)" << endl;
  cout << "Vega : " <<  call.vega() << " (LR)" << endl;
  cout << "Vega : " <<  call.vega("pw") << " (PW)" << endl;

    return 0;
}

int debug_konstantin(){
    european_option call(100, 100, 0.05, 0.4, 1, 0, 1000);
    cout<<"Theoretic price: "<<call.payoff_theoretic()<<endl;
    cout<<"Theoretic delta: "<<call.delta("th")<<endl;
    cout<<"Theoretic gamma: "<<call.gamma("th")<<endl;
    cout<<" Theoretic vega: "<<call.vega("th")<<endl;

    cout << "Simulation:"<<endl;
    cout << "Price: " << call.price() << endl;
    cout << "Delta: " << call.delta() << " (LR)" << endl;
    cout << "Delta: " << call.delta("pw") << " (PW)" << endl;
    cout << "Gamma: " << call.gamma() << " (PWLR)" << endl;
    cout << "Gamma: " << call.gamma("lrpw") << " (LRPW)" << endl;
    cout << "Gamma: " << call.gamma("lrlr") << " (LRLR)" << endl;
    cout << "Vega : " <<  call.vega() << " (LR)" << endl;
    cout << "Vega : " <<  call.vega("pw") << " (PW)" << endl;

    return 0;
}

int main(){
    cout << "hello" << endl;
    debug_thomas(1000);
    debug_thomas(10000);
    debug_thomas(100000);
    return 0;
}
