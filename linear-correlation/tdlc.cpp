#include<cmath>
#include<iostream>
#include<fstream>
#include<numeric>
#include<cstdlib>
#include"tdlc.h"

using namespace std;

void ReadData(string filename, vector<double> & data) {
	data.clear();
	ifstream ifile;
	ifile.open(filename.c_str());
	string s;
	string::size_type pos;
	string ss;
	while (getline(ifile, s)) {
		pos = s.find_first_of('\n', 0);
		ss = s.substr(0, pos);
		data.push_back(atof(ss.c_str()));
	}
	ifile.close();
}

double Mean(vector<int>& x) {
  double sum = accumulate(x.begin(), x.end(), 0.0);
  return sum * 1.0 / x.size();
}

double Mean(vector<double>& x) {
  double sum = accumulate(x.begin(), x.end(), 0.0);
  return sum / x.size();
}

double Std(vector<int>& x) {
  double x_mean = Mean(x);
  double std = 0;
  for (vector<int>::iterator it = x.begin(); it != x.end(); it ++) {
    std += pow(*it - x_mean, 2);
  }
  return sqrt(std);
}

double Std(vector<double>& x) {
  double x_mean = Mean(x);
  double std = 0;
  for (vector<double>::iterator it = x.begin(); it != x.end(); it ++) {
    std += pow(*it - x_mean, 2);
  }
  return sqrt(std);
}

double LC(vector<int>& raster, vector<double>& lfp) {
  double raster_mean = Mean(raster);
  double lfp_mean = Mean(lfp);
  vector<double> products(raster.size());
  for (int i = 0; i < raster.size(); i ++) {
    products[i] = (raster[i] - raster_mean) * (lfp[i] - lfp_mean);
  }
  return Mean(products) / (Std(raster) * Std(lfp));
}

void TDLC(vector<int>& raster, vector<double>& lfp, int negative_time_delay, int positive_time_delay, vector<double>& tdlc) {
  tdlc.resize(positive_time_delay + negative_time_delay + 1);
  vector<int> raster_copy = raster;
  vector<double> lfp_copy = lfp;
  for (int i = 0; i < negative_time_delay; i++) {
    raster_copy.erase(raster_copy.begin());
    lfp_copy.erase(lfp_copy.end() - 1);
    tdlc[i] = LC(raster_copy, lfp_copy);
  }
  tdlc[negative_time_delay] = LC(raster, lfp);
  raster_copy = raster;
  lfp_copy = lfp;
  for (int i = negative_time_delay + 1; i < negative_time_delay + positive_time_delay + 1; i++) {
    raster_copy.erase(raster_copy.end() - 1);
    lfp_copy.erase(lfp_copy.begin());
    tdlc[i] = LC(raster_copy, lfp_copy);
  }
}
