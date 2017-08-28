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
	//	STRING path: path of target file;
	//	VECTOR<VECTOR<T> > data: container of data;
	//	Return: none;
template <class T> void Read2D(string path, vector<vector<T> >& data) {
	data.clear();
	ifstream ifile;
	ifile.open(path.c_str());
	string s;
	vector<T> add_T;
	string::size_type last_pos = 0, current_pos = 0;
	while (getline(ifile, s)) {
		add_T.clear();
		current_pos = s.find_first_of(',', last_pos);
		while (current_pos != s.npos) {
			if (s.find('.') == string::npos) add_T.push_back(atoi(s.c_str()) + last_pos);
			else add_T.push_back(atof(s.c_str()) + last_pos);
			last_pos = current_pos + 1;
			current_pos = s.find_first_of(',', last_pos);
		}
		current_pos = s.find_first_of('\n', last_pos);
		if (current_pos != last_pos) {
			if (s.find('.') == string::npos) add_T.push_back(atoi(s.c_str()) + last_pos);
			else add_T.push_back(atof(s.c_str()) + last_pos);
		}
		data.push_back(add_T);
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
  	string s;
    string::size_type last_pos = 0, current_pos = 0;
  	int getline_counter = 0;
  	while (getline(ifile, s)) {
  		if (getline_counter == index) {
  			current_pos = s.find_first_of(',', last_pos);
  			while (current_pos != s.npos) {
  				if (s.find('.') == string::npos) data.push_back(atoi(s.c_str()) + last_pos);
					else data.push_back(atof(s.c_str()) + last_pos);
					last_pos = current_pos + 1;
  				current_pos = s.find_first_of(',', last_pos);
  			}
  			current_pos = s.find_first_of('\n', last_pos);
  			if (current_pos != last_pos) {
					if (s.find('.') == string::npos) data.push_back(atoi(s.c_str()) + last_pos);
					else data.push_back(atof(s.c_str()) + last_pos);
  			}
  			break;
  		}
  		getline_counter ++;
  	}
  } else {
    // prepare input file stream;
  	string s;
		string::size_type last_pos = 0, current_pos = 0;
  	int column_counter;
    bool whether_out;
  	while (getline(ifile, s)) {
      current_pos = s.find_first_of(',', last_pos);
      whether_out = false;
			column_counter = 0;
			while (current_pos != s.npos) {
  			if (column_counter == index) {
					if (s.find('.') == string::npos) data.push_back(atoi(s.c_str()) + last_pos);
					else data.push_back(atof(s.c_str()) + last_pos);
          whether_out = true;
  				break;
  			}
				last_pos = current_pos + 1;
        current_pos = s.find_first_of(',', last_pos);
        column_counter ++;
  		}
      if (whether_out == false) {
				if (column_counter == index) {
	        current_pos = s.find_first_of('\n', last_pos);
	    		if (current_pos != last_pos) {
						if (s.find('.') == string::npos) data.push_back(atoi(s.c_str()) + last_pos);
						else data.push_back(atof(s.c_str()) + last_pos);
	    		}
				}
      }
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
