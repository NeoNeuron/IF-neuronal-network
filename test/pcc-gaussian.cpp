#include<iostream>
#include<fstream>
#include<cstdlib>
#include<iomanip>
#include<cmath>
#include"../tdlc.h"
#include"../../tdmi/mi.h"
using namespace std;

int main(int argc, const char* argv[]) {
  if (argc != 2) throw runtime_error("wrong number of args");
	int length = atoi(argv[1]);
  vector<double> X(length), Y(length), U(length), V(length);
  for (int i = 0; i < length; i++){
    if (i == 0) {
      X[i] = GaussKernel();
      Y[i] = GaussKernel();
    } else {
      X[i] = -0.2 * Y[i - 1] + GaussKernel();
      Y[i] = -0.2 * X[i - 1] + GaussKernel();
    }
    U[i] = pow(X[i], 3);
    V[i] = Y[i];
  }
  // Main calculation:
  int negative_time_delay = 20;
  int positive_time_delay = 20;
  vector<double> tdlc_xy, tdlc_uv;
  // manipulate data;

  // cout << mean_first.size() << ',' << mean_second.size() << endl;
  TDLC(X, Y, negative_time_delay, positive_time_delay, tdlc_xy);
  TDLC(U, V, negative_time_delay, positive_time_delay, tdlc_uv);

  // output results:
  ofstream ofile;
  ofile.open("../file-txt/gaussian.txt");
  for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i ++) {
    ofile << (double) i - negative_time_delay << ',' << setprecision(6) << (double)tdlc_xy[i] << ',' << setprecision(6) << (double)tdlc_uv[i];
    if (i != (negative_time_delay + positive_time_delay)) ofile << endl;
  }
  ofile.close();
  return 0;
}
