//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-06-12
//	Description: Library for Data I/O functions in project;
//***************
#ifndef _IO_H_
#define _IO_H_

#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iomanip>
using namespace std;

//	Read 2 dimensional information;
	//	Read 2-D data from *.csv files; data type can be int or double;
	//	Return: none;
template <class T> void Read2D(string path, vector<vector<T> >& data) {
	data.clear();
	ifstream ifile;
	ifile.open(path.c_str());
	string s, ss;
	vector<T> add_T;
	string::size_type pos;
	while (getline(ifile, s)) {
		add_T.clear();
		pos = s.find_first_of(',', 0);
		while (pos != s.npos) {
			ss = s.substr(0, pos);
			if (ss.find('.') == string::npos) add_T.push_back(atoi(ss.c_str()));
			else add_T.push_back(atof(ss.c_str()));
			s.erase(0, pos + 1);
			ss.clear();
			pos = s.find_first_of(',', 0);
		}
		pos = s.find_first_of('\n', 0);
		if (pos == 0) continue;
		else {
			ss = s.substr(0, pos);
			if (ss.find('.') == string::npos) add_T.push_back(atoi(ss.c_str()));
			else add_T.push_back(atof(ss.c_str()));
		}
		data.push_back(add_T);
		s.clear();
	}
	ifile.close();
}

//	Read 1 dimensional (line/column) from *.csv file;
  //	STRING path: file path of *.csv file;
  //	INT index: the index of chosen line/column, from 0 to maximum - 1;
  //  INT axis: 0 for line, 1 for column;
  //	VECTOR<INT/DOUBLE> data: storage of output data;
  //	Return: none;
template <class T> void Read1D(string path, int index, int axis, vector<T>& data) {
  data.clear();
  // open data file;
  ifstream ifile;
  ifile.open(path.c_str());
  if (axis == 0) {
  	// prepare input file stream;
  	string s, ss;
    string::size_type pos;
    data.clear();
  	int getline_counter = 0;
  	while (getline(ifile, s)) {
  		if (getline_counter == index) {
  			pos = s.find_first_of(',', 0);
  			while (pos != s.npos) {
  				ss = s.substr(0, pos);
					if (ss.find('.') == string::npos) data.push_back(atoi(ss.c_str()));
					else data.push_back(atof(ss.c_str()));
  				s.erase(0, pos + 1);
  				ss.clear();
  				pos = s.find_first_of(',', 0);
  			}
  			pos = s.find_first_of('\n', 0);
  			if (pos != 0) {
  				ss = s.substr(0, pos);
					if (ss.find('.') == string::npos) data.push_back(atoi(ss.c_str()));
					else data.push_back(atof(ss.c_str()));
  			}
  			break;
  		}
      s.clear();
  		getline_counter ++;
  	}
  } else {
    // prepare input file stream;
  	string s, ss;
  	string::size_type pos;
  	int column_counter;
    bool whether_out;
  	while (getline(ifile, s)) {
      pos = s.find_first_of(',', 0);
      whether_out = false;
			column_counter = 0;
			while (pos != s.npos) {
  			if (column_counter == index) {
  				ss = s.substr(0, pos);
					if (ss.find('.') == string::npos) data.push_back(atoi(ss.c_str()));
					else data.push_back(atof(ss.c_str()));
          whether_out = true;
  				break;
  			}
  			s.erase(0, pos + 1);
        pos = s.find_first_of(',', 0);
        column_counter ++;
  		}
      if (whether_out == false) {
				if (column_counter == index) {
	        pos = s.find_first_of('\n', 0);
	    		if (pos == 0) continue;
	    		else {
	    			ss = s.substr(0, pos);
						if (ss.find('.') == string::npos) data.push_back(atoi(ss.c_str()));
						else data.push_back(atof(ss.c_str()));
	    		}
				}
      }
      s.clear();
  	}
  }
  ifile.close();
}

//  Print 2 dimensional data;
  //  Print 2-D data to *.csv files; data type can be int, double or bool;
  //  STRING mode: openmode for aiming file; "app" for append, "trunc" for any contents that existed in the file before it is open are discarded.
  //  Return: none;
template <class T> void Print2D(string path, string mode, vector<vector<T> >& data) {
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

//	Print 1 dimensional (line/column) to *.csv file;
  //	STRING path: file path of *.csv file;
  //  STRING mode: openmode for aiming file; "app" for append, "trunc" for any contents that existed in the file before it is open are discarded.
  //  INT axis: 0 for line, 1 for column;
  //	VECTOR<INT/DOUBLE> data: storage of outputing data;
  //	Return: none;
template <class T> void Print1D(string path, string mode, int axis, vector<T>& data) {
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
#endif // _IO_H_
