#include <fstream>
#include <ctime>

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
