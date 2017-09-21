//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-08-29
//	Description: Library for Data I/O functions in project;
//***************
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;

//	Read 2 dimensional information;
	//	Read 2-D data from *.csv files; data type can be int or double;
	//	STRING path: path of target file;
	//	VECTOR<VECTOR<T> > data: container of data;
	//	Return: none;
template <class T> void Read2D(string path, vector<vector<T> >& data) {
  data.clear();
	ifstream ifile(path.c_str());
	string s, buffer;
	vector<T> add_T;
	while (getline(ifile, s)) {
		add_T.clear();
    istringstream s_input(s);
		while (getline(s_input, buffer, ',')) {
      if (buffer.find('.') == string::npos) {
        add_T.push_back(atoi(buffer.c_str()));
      } else add_T.push_back(atof(buffer.c_str()));
		}
		data.push_back(add_T);
	}
	ifile.close();
}

template <class T> void Read1D(string path, vector<T>& data, int index, int axis) {
  data.clear();
  // open data file;
  ifstream ifile(path.c_str());
  if (axis == 0) {
  	string s, buffer;
  	int getline_counter = 0;
  	while (getline(ifile, s)) {
  		if (getline_counter == index) {
  			istringstream s_input(s);
  			while (getline(s_input, buffer, ',')) {
  				if (buffer.find('.') == string::npos) {
            data.push_back(atoi(buffer.c_str()));
					} else data.push_back(atof(buffer.c_str()));
  			}
  		}
  		getline_counter ++;
  	}
  } else {
  	string s, buffer;
  	int column_counter;
  	while (getline(ifile, s)) {
			column_counter = 0;
      istringstream s_input(s);
			while (getline(s_input, buffer, ',')) {
  			if (column_counter == index) {
					if (buffer.find('.') == string::npos) {
            data.push_back(atoi(buffer.c_str()));
					} else data.push_back(atof(buffer.c_str()));
  				break;
  			}
        column_counter ++;
  		}
  	}
  }
  ifile.close();
}

template <class T> void Print2D(string path, vector<vector<T> >& data, string mode) {
	ofstream ofile;
	if (mode == "app") ofile.open(path.c_str(), ios::app);
	else if (mode == "trunc") ofile.open(path.c_str());
	for (typename vector<vector<T> >::iterator it = data.begin(); it != data.end(); it++) {
		for (typename vector<T>::iterator itt = it->begin(); itt != it->end(); itt++) {
			ofile << setprecision(15) << *itt << ',';
		}
		ofile << '\n';
	}
	ofile.close();
}

template <class T> void Print1D(string path, vector<T>& data, string mode, int axis) {
	ofstream ofile;
	if (mode == "app") ofile.open(path.c_str(), ios::app);
	else if (mode == "trunc") ofile.open(path.c_str());
	if (axis == 0) {
		for (typename vector<T>::iterator it = data.begin(); it != data.end(); it++) {
			ofile << setprecision(15) << *it << ',';
		}
		ofile << '\n';
	} else {
		for (typename vector<T>::iterator it = data.begin(); it != data.end(); it++) {
			ofile << setprecision(15) << *it << '\n';
		}
	}
	ofile.close();
}
