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
#include <cmath>
#include <algorithm>

using namespace std;


int main(){
    cout << "hello" << endl;




    lookback_option l(160, 130, 0.05, 0.4, 1, 0, 1000000);
    cout << "lookback LR delta  : " << l.delta("lr") << endl;
    cout << "lookback PW delta  : " << l.delta("pw") << endl;
    cout << "lookback TH delta  : " << l.delta("th") << endl;
    cout << "lookback PWLR gamma: " << l.gamma("pwlr") << endl;
    cout << "lookback TH gamma  : " << l.gamma("th") << endl;
    cout << "lookback LR  vega  : " << l.vega("lr") << endl;
    cout << "lookback PW  vega  : " << l.vega("pw") << endl;
    cout << "lookback TH vega   : " << l.vega("th") << endl;
    cout << "lookback    Price  : " << l.price() << endl;
    cout << "lookback TH price  : " << l.payoff_theoretic() << endl;

    return 0;
}

