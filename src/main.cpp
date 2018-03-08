/*
This is the main project file.
*/

#include<iostream> // For console output
#include <fstream> // For csv output

#include "european_option.hpp"
#include "barrier_option.hpp"


#include <cmath> // For math calculations
#include <string> // For user instructions in greeks
#include <random> // For random numbers generation

using namespace std;


int write_csv(int pace){
  ofstream myfile;
  myfile.open("../report/graphs/core_task.csv");
  myfile << "iterations;delta lr;delta pw;gamma pwlr;gamma lrpw;gamma lrlr;vega lr;vega pw\n";

  european_option s1(100, 100, 0.05, 0.4, 1, 0, 1);

  for (int i=1; i<10000000; i*=pace){
    s1 = european_option(100, 100, 0.05, 0.4, 1, 0, i);
    myfile << i << ";" << s1.delta() << ";" << s1.delta("lr") << ";" << s1.gamma() << ";" << s1.gamma("lrpw") << ";" << s1.gamma("lrlr") << ";" << s1.vega() << ";" << s1.vega("pw") << "\n";
  }


  return 0;
}



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
    //write_csv(2);
    barrier_option boption(90, 100, 100, 0.05, 0.4, 1, 0, 1000);
    cout << "Price: " << boption.price() << endl;
    cout << "Price: " << boption.price_theoretic() << " (Theoretical)" << endl;
    cout << "Delta: " << boption.delta() << " (LR)" << endl;
    cout << "Delta: " << boption.delta("th") << " (Theoretical)" << endl;
    cout << "Gamma: " << boption.gamma() << " (LR)" << endl;
    cout << "Gamma: " << boption.gamma("th") << " (Theoretical)" << endl;
    cout << "Vega: " << boption.vega() << " (LR)" << endl;
    cout << "Vega: " << boption.vega("th") << " (Theoretical)" << endl;
    return 0;
}
