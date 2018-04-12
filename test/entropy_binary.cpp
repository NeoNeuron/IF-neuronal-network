#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cmath>
#include <algorithm>
#include "../include/io.h"
using namespace std;

int Max(vector<int>& data) {
	return *(max_element(data.begin(), data.end()));
}

void P(vector<int>& x, vector<int>& px) {
	// find the maximum value of x;
	int xmax = Max(x);
	px.clear();
	px.resize(xmax+1, 0);
	for (unsigned int i = 0; i < x.size(); i++) {
		px[x[i]] ++;
	}
	// cout << "1-2" << endl;
	return;
}

double H(vector<int>& x) {
	vector<int> px(x.size());
	P(x, px);
	double H_var = 0.0;
	for (int i = 0; i < px.size(); i++) {
		if (px[i] != 0) H_var += px[i]*log2(px[i]);
	}
	// cout << "1-3" << endl;
	return -H_var/x.size() + log2(x.size());
}

void bin2int(vector<bool>& bdata, vector<int>& idata, int k) {
	int N = floor(bdata.size()/k);
	idata.resize(N);
	long int tmp;
	for (int i = 0; i < N; i ++) {
		tmp = 0;
		for (int j = 0; j < k; j ++) {
			if (bdata[i*k + j]) {
				tmp += pow(2,j);
			}
		}
		idata[i] = tmp;
	}
	// cout << "1-1" << endl;
}

int main() {
	// Setups for experiment;
	// length of data;
	int length = 2e4;
	// model firing rate;
	double p = 0.02;
	// range of k:
	int kmax = 30;

	// preparing the model spike train;
	vector<bool> s(length);
	double xrand;
	for (int i = 0; i < length; i ++) {
		xrand = rand()/(RAND_MAX+1.0);
		if (xrand < p) s[i]= true;
		else s[i] = false;
	}
	// calculate the entropy of x;
	vector<double> hs(kmax);
	vector<int> idata;
	for (int i = 0; i < kmax; i ++) {
		cout << i + 1 << endl;
		idata.clear();
		bin2int(s, idata, i + 1);
		hs[i] = H(idata);
	}
	Print1D("./data/tmp/entro_bin.csv", hs, "trunc", 1);
	return 0;
}
