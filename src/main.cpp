/*
This is the main project file.
*/

#include "european_option.hpp"
#include "barrier_option.hpp"
#include "lookback_option.hpp"

#include <iostream> // For console output
#include <cmath> // For math calculations
#include <string> // For user instructions in greeks
#include <random> // For random numbers generation


#include <fstream>
#include <ctime>

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
    lookback_option l(3.0, 100, 120, 0.05, 0.3, 0, 100000, 100);
    cout << "lookback lR delta: " << l.delta() << endl;
    cout << "lookback PW delta: " << l.delta("pw") << endl;
    cout << "lookback LR gamma: " << l.gamma() << endl;
    cout << "lookback LR  vega: " << l.vega() << endl;
    cout << "lookback    Price: " << l.price() << endl;

    return 0;
}


int write_csv(string pathToFile){
  ofstream myfile;
  myfile.open(pathToFile);
  myfile << "iterations;delta th;delta th time;delta lr;delta lr time;delta pw;delta pw time;gamma th;gamma th time;gamma pwlr;gamma pwlr time;gamma lrpw;gamma lrpw time;gamma lrlr;gamma lrlr time;vega th;vega th time;vega lr;vega lr time;vega pw;vega pw time\n";

  european_option s1(100, 100, 0.05, 0.4, 1, 0, 1);

  for (int i=1; i<10000000; i*=10){
    for (int j=0; j<1000; j+=1){
      s1 = european_option(100, 100, 0.05, 0.4, 1, 0, i);

      clock_t begin = clock();
      double delta_th = s1.delta("th");
      clock_t end = clock();
      double delta_th_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double delta_lr = s1.delta("lr");
      end = clock();
      double delta_lr_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double delta_pw = s1.delta("pw");
      end = clock();
      double delta_pw_time = double(end - begin) / CLOCKS_PER_SEC;


      myfile << i << ";" << delta_th << ";" << delta_th_time;
      myfile << ";" << delta_lr << ";" << delta_lr_time;
      myfile << ";" << delta_pw << ";" << delta_pw_time;

      begin = clock();
      double gamma_th = s1.gamma("th");
      end = clock();
      double gamma_th_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double gamma_pwlr = s1.gamma("pwlr");
      end = clock();
      double gamma_pwlr_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double gamma_lrpw = s1.gamma("lrpw");
      end = clock();
      double gamma_lrpw_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double gamma_lrlr = s1.gamma("lrlr");
      end = clock();
      double gamma_lrlr_time = double(end - begin) / CLOCKS_PER_SEC;


      myfile << ";" << gamma_th << ";" << gamma_th_time;
      myfile << ";" << gamma_pwlr << ";" << gamma_pwlr_time;
      myfile << ";" << gamma_lrpw << ";" << gamma_lrpw_time;
      myfile << ";" << gamma_lrlr << ";" << gamma_lrlr_time;

      begin = clock();
      double vega_th = s1.vega("th");
      end = clock();
      double vega_th_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double vega_lr = s1.vega("lr");
      end = clock();
      double vega_lr_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double vega_pw = s1.vega("pw");
      end = clock();
      double vega_pw_time = double(end - begin) / CLOCKS_PER_SEC;

      myfile << ";" << vega_th << ";" << vega_th_time;
      myfile << ";" << vega_lr << ";" << vega_lr_time;
      myfile << ";" << vega_pw << ";" << vega_pw_time;


      myfile << "\n";

    }
  }


  return 0;
}



int write_csvlookback(string pathToFile){
  ofstream myfile;
  myfile.open(pathToFile);
  myfile << "iterations;delta lr;delta lr time;delta pw;delta pw time;gamma pwlr;gamma pwlr time;vega lr;vega lr time;\n";

  lookback_option s1(100, 100, 0.05, 0.4, 1, 0, 1);

  for (int i=1; i<10000000; i*=10){
    for (int j=0; j<1000; j+=1){
      s1 = lookback_option(100, 100, 0.05, 0.4, 1, 0, i);

   /*   clock_t begin = clock();
      double delta_th = s1.delta("th");
      clock_t end = clock();
      double delta_th_time = double(end - begin) / CLOCKS_PER_SEC;
     */
      clock_t begin = clock();
      double delta_lr = s1.delta("lr");
      clock_t end = clock();
      double delta_lr_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double delta_pw = s1.delta("pw");
      end = clock();
      double delta_pw_time = double(end - begin) / CLOCKS_PER_SEC;


      //myfile << i << ";" << delta_th << ";" << delta_th_time;
      myfile << ";" << delta_lr << ";" << delta_lr_time;
      myfile << ";" << delta_pw << ";" << delta_pw_time;

     /* begin = clock();
      double gamma_th = s1.gamma("th");
      end = clock();
      double gamma_th_time = double(end - begin) / CLOCKS_PER_SEC;
    */begin = clock();
      double gamma_pwlr = s1.gamma();
      end = clock();
      double gamma_pwlr_time = double(end - begin) / CLOCKS_PER_SEC;
     /* begin = clock();
      double gamma_lrpw = s1.gamma("lrpw");
      end = clock();
      double gamma_lrpw_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double gamma_lrlr = s1.gamma("lrlr");
      end = clock();
      double gamma_lrlr_time = double(end - begin) / CLOCKS_PER_SEC;
    */

      //myfile << ";" << gamma_th << ";" << gamma_th_time;
      myfile << ";" << gamma_pwlr << ";" << gamma_pwlr_time;
      //myfile << ";" << gamma_lrpw << ";" << gamma_lrpw_time;
      //myfile << ";" << gamma_lrlr << ";" << gamma_lrlr_time;

      begin = clock();
      double vega_th = s1.vega();
      end = clock();
      double vega_th_time = double(end - begin) / CLOCKS_PER_SEC;
      /*begin = clock();
      double vega_lr = s1.vega("lr");
      end = clock();
      double vega_lr_time = double(end - begin) / CLOCKS_PER_SEC;
      begin = clock();
      double vega_pw = s1.vega("pw");
      end = clock();
      double vega_pw_time = double(end - begin) / CLOCKS_PER_SEC;
*/
      myfile << ";" << vega_th << ";" << vega_th_time;
     // myfile << ";" << vega_lr << ";" << vega_lr_time;
      //myfile << ";" << vega_pw << ";" << vega_pw_time;


      myfile << "\n";

    }
  }


  return 0;
}

int main(){
    cout << "hello" << endl;
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
    //debug_konstantin();
    return 0;
}
