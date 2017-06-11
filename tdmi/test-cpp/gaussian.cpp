//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-03
//	Description: Apply Gaussian series to test the mutual information analysis
//***************
#include "../mi.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <cstdlib>

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
      X[i] = -0.2 * X[i - 1] + GaussKernel();
      Y[i] = -0.2 * X[i - 1] + GaussKernel();
    }
    U[i] = pow(X[i], 3);
    V[i] = Y[i];
  }
  // cout << "part1" << endl;
  // TDMI with a uniform partition of XY plane
  //	Calculate time-delayed mutual information;
  vector<double> xy_tdmi, uv_tdmi;
  int positive_time_delay = 20;
  int negative_time_delay = 20;
  TDMI_uniform(X, Y, negative_time_delay, positive_time_delay, xy_tdmi);
  TDMI_uniform(U, V, negative_time_delay, positive_time_delay, uv_tdmi);
  // Output data;
  ofstream tdmi_test;
  tdmi_test.open("../file-txt/gaussian_uniform.txt");
  for (int i = 0; i < positive_time_delay + negative_time_delay + 1; i++) {
    tdmi_test << i - negative_time_delay << '\t' << xy_tdmi[i] << "\t" << uv_tdmi[i];
    if (i < positive_time_delay + negative_time_delay) tdmi_test << endl;
  }
  tdmi_test.close();
  // cout << "part1" << endl;
  // TDMI with an adaptive partition of XY plane
	//	Calculate time-delayed mutual information;
	xy_tdmi.clear();
  uv_tdmi.clear();
	positive_time_delay = 20;
	negative_time_delay = 20;
	TDMI_adaptive(X, Y, negative_time_delay, positive_time_delay, xy_tdmi);
	TDMI_adaptive(U, V, negative_time_delay, positive_time_delay, uv_tdmi);

	tdmi_test.open("../file-txt/gaussian_adaptive.txt");
	for (int i = 0; i < positive_time_delay + negative_time_delay + 1; i++) {
		tdmi_test << i - negative_time_delay << '\t' << xy_tdmi[i] << "\t" << uv_tdmi[i];
		if (i < positive_time_delay + negative_time_delay) tdmi_test << endl;
	}
	tdmi_test.close();
	return 0;
}
