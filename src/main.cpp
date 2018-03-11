/*
This is the main project file.
*/
#include "barrierlookbackoptions.h"

#include<iostream> // For console output
using namespace std;

#include <ctime>




int debug_thomas(){
    barrier_option call1(80, 100, 100, 0.05, 0.4, 1, 0, 10000);
    lookback_option call2(80, 100, 100, 0.05, 0.4, 1, 0, 10000);
    //cout<<"Theoretic price: "<<call.payoff_theoretic()<<endl;
    //cout<<"Theoretic delta: "<<call.delta()<<endl;
    //cout<<"Theoretic gamma: "<<call.gamma()<<endl;


    clock_t begin = clock();
    call1.price();
    clock_t end = clock();
    cout << "BarrierOption Time: " << double(end - begin) / CLOCKS_PER_SEC << endl;

    begin = clock();
    call2.price();
    end = clock();
    cout << "LookBackOption Time: " << double(end - begin) / CLOCKS_PER_SEC << endl;
    /*
    cout << "Simulation:"<<endl;
    cout << "Price: " << call.price() << endl;
    cout << "Delta: " << call.delta() << " (LR)" << endl;
    cout << "Delta: " << call.delta("pw") << " (PW)" << endl;
    cout << "Gamma: " << call.gamma() << " (PWLR)" << endl;
    cout << "Gamma: " << call.gamma("lrpw") << " (LRPW)" << endl;
    cout << "Gamma: " << call.gamma("lrlr") << " (LRLR)" << endl;
    */
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

        barrier_option call(125, 100, 100, 0.05, 0.4, 1, 0, 10000);   // need Price  around 18
    cout << "Theoretic price: " << call.price_theoretic() << endl;
    cout << "Price: " << call.price() << endl;
    cout << "Delta: " << call.delta("lr") << " (LR)" << endl;
    cout << "Delta: " << call.delta("th") << " (th)" << endl;
    cout << "Gamma: " << call.gamma("lr") << " (LRLR)" << endl;
    cout << "Gamma: " << call.gamma("th") << " (th)" << endl;
    cout << "Vega : " <<  call.vega("lr") << " (LR)" << endl;
    cout << "Vega : " <<  call.vega("th") << " (th)" << endl;

*/

    lookback_option l(100, 120, 0.05, 0.3, 3.0, 0, 100000, 100);
    cout << "lookback LR delta  : " << l.delta("lr") << endl;
    cout << "lookback PW delta  : " << l.delta("pw") << endl;
    cout << "lookback TH delta  : " << l.delta("th") << endl;
    cout << "lookback LRLR gamma: " << l.gamma("lrlr") << endl;
    cout << "lookback PWLR gamma: " << l.gamma("pwlr") << endl;
    cout << "lookback TH gamma  : " << l.gamma("th") << endl;
    cout << "lookback LR  vega  : " << l.vega("lr") << endl;
    cout << "lookback PW  vega  : " << l.vega("pw") << endl;
    cout << "lookback TH vega   : " << l.vega("th") << endl;
    cout << "lookback    Price  : " << l.price() << endl;
    cout << "lookback TH price  : " << l.payoff_theoretic() << endl;

    return 0;
}

int main(){
    cout << "hello" << endl;
    //debug_konstantin();
    debug_thomas();
    return 0;
}
