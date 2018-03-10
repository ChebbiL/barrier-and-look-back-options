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
    /*
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
    lookback_option l(100, 120, 0.05, 0.3, 3.0, 0, 100000, 100);
    cout << "lookback lR delta: " << l.delta() << endl;
    cout << "lookback PW delta: " << l.delta("pw") << endl;
    cout << "lookback LR gamma: " << l.gamma() << endl;
    cout << "lookback LR  vega: " << l.vega() << endl;
    cout << "lookback    Price: " << l.price() << endl;
*/
    barrier_option call(80, 100, 100, 0.05, 0.4, 1, 0, 10000);   // need Price  around 18
    cout << "Theoretic price: " << call.price_theoretic() << endl;
    cout << "Price: " << call.price() << endl;
    cout << "Delta: " << call.delta("lr") << " (LR)" << endl;
    cout << "Delta: " << call.delta("th") << " (th)" << endl;
    cout << "Gamma: " << call.gamma("lr") << " (LRLR)" << endl;
    cout << "Gamma: " << call.gamma("th") << " (th)" << endl;
    cout << "Vega : " <<  call.vega("lr") << " (LR)" << endl;
    cout << "Vega : " <<  call.vega("th") << " (th)" << endl;

    return 0;
}

int main(){
    cout << "hello" << endl;
<<<<<<< HEAD
    /*barrier_option boption(50, 100, 100, 0.05, 0.4, 1, 0, 1000);
    cout << "Price: " << boption.price() << endl;
    cout << "Price: " << boption.price_theoretic() << " (Theoretical)" << endl;
    cout << "Delta: " << boption.delta() << " (LR)" << endl;
    cout << "Delta: " << boption.delta("th") << " (Theoretical)" << endl;
    cout << "Gamma: " << boption.gamma() << " (LR)" << endl;
    cout << "Gamma: " << boption.gamma("th") << " (Theoretical)" << endl;
    cout << "Vega: " << boption.vega() << " (LR)" << endl;
    cout << "Vega: " << boption.vega("th") << " (Theoretical)" << endl;
*/
    //debug_thomas();
    //write_csv("../report/graphs/core_task.csv");
    write_csvlookback("../report/graphs/lookback.csv");
    //debug_konstantin();
=======
    debug_konstantin();

>>>>>>> c2c060cc2739de0cbbb9d57fec211a40872e76fc
    return 0;
}
