//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-02-27 21:47:53
//	Description: Test all four TDMI functions;
//***************
#include "../mi.h"
#include "../../io/io.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <ctime>

using namespace std;

int SetOption() {
	cout << "Choose the function you want to test: " << endl;
	cout << ">> [1]: TDMI between two spike trains;" << endl;
	cout << ">> [2]: TDMI between two continues double sequence; applying uniformly-binned historgram scheme;" << endl;
	cout << ">> [3]: TDMI between two continues double sequence, applying histogram scheme with equal occupancy for each bin;" << endl;
	cout << ">> [4]: TDMI between spike train and local field potential, applying uniformly-binned historgram scheme;" << endl;
	cout << ">> [5]: TDMI between spike train and local field potential, applying histogram scheme with equal occupancy for each bin;" << endl;
	cout << ">> [0]: Exit;" << endl;
	int option;
	cin >> option;
	return option;
}

int main() {
	int option;
	option = SetOption();
	while (option != 0) {
		if (option == 1) {
			double tmax = 10000;
			double x_rate = 5, y_rate = 5; // Poisson rate of x and y;
			vector<double> x, y;
			x.push_back(0);
			y.push_back(0);
			double new_x = 0, new_y = 0;
			while (new_x < tmax) {
				new_x = PoissonGenerator(x_rate, x.back());
				x.push_back(new_x);
			}
			while (new_y < tmax) {
				new_y = PoissonGenerator(y_rate, y.back());
				y.push_back(new_y);
			}
			//	Calculate time-delayed mutual information;
			vector<double> xy_tdmi;
			int positive_time_delay = 20;
			int negative_time_delay = 20;
			double dt = 0.125;
			TDMI(x, y, dt, tmax, 1, negative_time_delay, positive_time_delay, xy_tdmi);
			// Output data;
			ofstream tdmi_test;
			tdmi_test.open("tdmi_test_1.dat");
			for (int i = 0; i < positive_time_delay + negative_time_delay + 1; i++) {
				tdmi_test << (double)dt*(i - negative_time_delay) << ',' << xy_tdmi[i];
				if (i < positive_time_delay + negative_time_delay) tdmi_test << endl;
			}
			tdmi_test.close();
			option = SetOption();
		} else if (option == 2) {
			int set_size = 10000; // size of x, y, u, v;
			vector<double> x(set_size, 0), y(set_size, 0), u(set_size, 0), v(set_size, 0); // x, y are two int random variables with Gaussion distribution; u = f(x), v = g(y);
			x[0] = GaussKernel();
			y[0] = GaussKernel();
			double correlation_factor = 0.1;
			for (int i = 1; i < set_size; i++) {
				x[i] = -correlation_factor*y[i-1] + GaussKernel();
				y[i] = -correlation_factor*x[i-1] + GaussKernel();
			}
			for (int i = 0; i < set_size; i++) {
				u[i] = (x[i]);
				v[i] = (abs(y[i]) * y[i]);
			}

			//	Calculate time-delayed mutual information;
			vector<double> xy_tdmi, uv_tdmi;
			int positive_time_delay = 20;
			int negative_time_delay = 20;
			TDMI(x, y, 0.1, 0.1, negative_time_delay, positive_time_delay, xy_tdmi);
			TDMI(u, v, 0.1, 0.1, negative_time_delay, positive_time_delay, uv_tdmi);
			// Output data;
			ofstream tdmi_test;
			tdmi_test.open("tdmi_test_2.dat");
			for (int i = 0; i < positive_time_delay + negative_time_delay + 1; i++) {
				tdmi_test << i - negative_time_delay << ',' << xy_tdmi[i] << "," << uv_tdmi[i];
				if (i < positive_time_delay + negative_time_delay) tdmi_test << endl;
			}
			tdmi_test.close();
			option = SetOption();
		} else if (option == 3) {
			int set_size = 10000; // size of x, y, u, v;
			vector<double> x(set_size, 0), y(set_size, 0), u(set_size, 0), v(set_size, 0); // x, y are two int random variables with Gaussion distribution; u = f(x), v = g(y);
			x[0] = GaussKernel();
			y[0] = GaussKernel();
			double correlation_factor = 0.1;
			for (int i = 1; i < set_size; i++) {
				x[i] = -correlation_factor*y[i-1] + GaussKernel();
				y[i] = -correlation_factor*x[i-1] + GaussKernel();
			}
			for (int i = 0; i < set_size; i++) {
				u[i] = (x[i]);
				v[i] = (abs(y[i]) * y[i]);
			}

			//	Calculate time-delayed mutual information;
			int expected_occupancy = 10;
			vector<double> xy_tdmi, uv_tdmi;
			int positive_time_delay = 20;
			int negative_time_delay = 20;
			TDMI(x, y, expected_occupancy, negative_time_delay, positive_time_delay, xy_tdmi);
			TDMI(u, v, expected_occupancy, negative_time_delay, positive_time_delay, uv_tdmi);

			ofstream tdmi_test;
			tdmi_test.open("tdmi_test_3.dat");
			for (int i = 0; i < positive_time_delay + negative_time_delay + 1; i++) {
				tdmi_test << i - negative_time_delay << ',' << xy_tdmi[i] << "," << uv_tdmi[i];
				if (i < positive_time_delay + negative_time_delay) tdmi_test << endl;
			}
			tdmi_test.close();
			option = SetOption();
		} else if (option == 4) {
			// INPUT NEURONAL DATA:
			vector<double> raster_test;
			vector<double> LFP_test;
			string folder = "./test_file/";
			ifstream data;
			string path;
			path = folder + "LFP_test.txt";
			Read1D(path, 0, 1, LFP_test);
			path = folder + "raster_test.txt";
			Read1D(path, 0, 1, raster_test);

			// TIME DELAY
			vector<double> tdmi;
			double dt = 0.125;
			int positive_time_delay = 20;
			int negative_time_delay = 20;
			TDMI(raster_test, LFP_test, dt, 0.03125, 0.01, negative_time_delay, positive_time_delay, tdmi);

			ofstream tdmi_test;
			tdmi_test.open("tdmi_test_4.dat");
			for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
				tdmi_test << (double) dt*(i - negative_time_delay) << ',' << (double)tdmi[i];
				if (i < positive_time_delay + negative_time_delay) tdmi_test << endl;
			}
			tdmi_test.close();
			option = SetOption();
		} else if (option == 5) {
			// INPUT NEURONAL DATA:
			vector<double> raster_test;
			vector<double> LFP_test;
			string folder = "./test_file/";
			ifstream data;
			string path;
			path = folder + "LFP_test.txt";
			Read1D(path, LFP_test);
			path = folder + "raster_test.txt";
			Read1D(path, raster_test);

			// TIME DELAY
			vector<double> tdmi;
			int expected_occupancy = 25;
			double dt = 0.125;
			int positive_time_delay = 80;
			int negative_time_delay = 100;
			TDMI(raster_test, LFP_test, expected_occupancy, dt, 0.03125, negative_time_delay, positive_time_delay, tdmi, false);

			ofstream tdmi_test;
			tdmi_test.open("tdmi_test_5.dat");
			for (int i = 0; i < negative_time_delay + positive_time_delay + 1; i++) {
				tdmi_test << (double)dt*(i - negative_time_delay) << ',' << (double)tdmi[i];
				if (i < positive_time_delay + negative_time_delay) tdmi_test << endl;
			}
			tdmi_test.close();
			option = SetOption();
		} else {
			option = SetOption();
			continue;
		}
	}
	return 0;
}
