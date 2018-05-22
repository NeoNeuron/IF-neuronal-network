//*************************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-04-03 15:07:52
//	Description: test program for multi-network simulation;
//*************************

#include "../include/connectivity_matrix.h"
#include "../include/io.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>
using namespace std;

//	Calculate the least path matrix for given network structure;
//	arguments:
//	argv[1] = path of connectivity matrix;
//	argv[2] = path of outputing least path matrix;
int main(int argc, const char* argv[]) {
	if (argc != 3) throw runtime_error("wrong number of args");
	clock_t start, finish;
	start = clock();
	// 	Setup directory for output files;
	//	it must be existing dir;
  string mat_dir = argv[1];
  string path_mat_dir = argv[2];

	// Loading matfile:
	vector<vector<int> > mat_data;
  Read2D(mat_dir, mat_data);

	// load neuron number;
	int neuron_number = 100;
	ConnectivityMatrix mat;
  mat.SetNeuronNumber(neuron_number);

  mat.LoadMatrix(mat_data);
  mat.FindLeastPath();
  mat.OutPathMatrix(path_mat_dir);

	finish = clock();

	// COUNTS:
	double ToTtime;
	ToTtime = (finish - start)*1.0 / CLOCKS_PER_SEC;
	printf(">> It takes %.2fs\n", ToTtime);
	return 0;
}
