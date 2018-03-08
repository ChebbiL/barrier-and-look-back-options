/*
This is the main project file.
*/

#include "european_option.hpp"
#include "barrier_option.hpp"
#include "lookback_option.hpp"

#include<iostream> // For console output
#include <cmath> // For math calculations
#include <string> // For user instructions in greeks
#include <random> // For random numbers generation

using namespace std;





int debug_thomas(){
    european_option call(100, 100, 0.05, 0.4, 1, 0, 1000);
    cout<<"Theoretic price: "<<call.payoff_theoretic()<<endl;
    cout<<"Theoretic delta: "<<call.delta("th")<<endl;
    cout<<"Theoretic gamma: "<<call.gamma("th")<<endl;

    cout << "Simulation:"<<endl;
    cout << "Price: " << call.price() << endl;
    cout << "Delta: " << call.delta() << " (LR)" << endl;
    cout << "Delta: " << call.delta("pw") << " (PW)" << endl;
    cout << "Gamma: " << call.gamma() << " (PWLR)" << endl;
    cout << "Gamma: " << call.gamma("lrpw") << " (LRPW)" << endl;
    cout << "Gamma: " << call.gamma("lrlr") << " (LRLR)" << endl;

    return 0;
}

int debug_konstantin(){
    cout<<fixed;
    european_option call(100, 100, 0.05, 0.4, 1, 0, 1000000);
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
    cout<< "______________________________"<<endl;

    lookback_option l(100, 120, 0.05, 0.3, 3.0, 0, 100000, 100, true);
    cout << "lookback lR delta: " << l.delta_LR() << endl;
    cout << "lookback PW delta: " << l.delta_PW() << endl;
    cout << "lookback LR gamma: " << l.gamma_LR() << endl;
    cout << "lookback LR  vega: " << l.vega_LR() << endl;
    cout << "lookback    Price: " << l.value() << endl;

    return 0;
}

int main(){
    cout << "hello" << endl;
    barrier_option boption(90, 100, 100, 0.05, 0.4, 1, 0, 1000);
    cout << "Price: " << boption.price() << endl;
    cout << "Price: " << boption.payoff_theoretic() << " (Theoretical)" << endl;
    cout << "Delta: " << boption.delta() << " (LR)" << endl;
    cout << "Delta: " << boption.delta("th") << " (Theoretical)" << endl;
    cout << "Gamma: " << boption.gamma() << " (LR)" << endl;
    cout << "Gamma: " << boption.gamma("th") << " (Theoretical)" << endl;
    cout << "Vega: " << boption.vega() << " (LR)" << endl;
    cout << "Vega: " << boption.vega("th") << " (Theoretical)" << endl;
    //write_csv(2);
    debug_konstantin();
    return 0;
}

