#include "io.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>

using namespace std;

void Read2D(string path, vector<vector<int> >& data) {
	data.clear();
	ifstream ifile;
	ifile.open(path.c_str());
	string s, ss;
	vector<int> add_int;
	string::size_type pos;
	while (getline(ifile, s)) {
		add_int.clear();
		pos = s.find_first_of(',', 0);
		while (pos != s.npos) {
			ss = s.substr(0, pos);
			add_int.push_back(atoi(ss.c_str()));
			s.erase(0, pos + 1);
			ss.clear();
			pos = s.find_first_of(',', 0);
		}
		pos = s.find_first_of('\n', 0);
		if (pos == 0) continue;
		else {
			ss = s.substr(0, pos);
			add_int.push_back(atoi(ss.c_str()));
		}
		data.push_back(add_int);
		s.clear();
	}
	ifile.close();
}

void Read2D(string path, vector<vector<double> >& data) {
	data.clear();
	ifstream ifile;
	ifile.open(path.c_str());
	string s, ss;
	vector<double> add_double;
	string::size_type pos;
	while (getline(ifile, s)) {
		add_double.clear();
		pos = s.find_first_of(',', 0);
		while (pos != s.npos) {
			ss = s.substr(0, pos);
			add_double.push_back(atof(ss.c_str()));
			s.erase(0, pos + 1);
			ss.clear();
			pos = s.find_first_of(',', 0);
		}
		pos = s.find_first_of('\n', 0);
		if (pos == 0) continue;
		else {
			ss = s.substr(0, pos);
			add_double.push_back(atof(ss.c_str()));
		}
		data.push_back(add_double);
		s.clear();
	}
	ifile.close();
}

void Read1D(string path, int index, int axis, vector<int>& data) {
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
  				data.push_back(atoi(ss.c_str()));
  				s.erase(0, pos + 1);
  				ss.clear();
  				pos = s.find_first_of(',', 0);
  			}
  			pos = s.find_first_of('\n', 0);
  			if (pos != 0) {
  				ss = s.substr(0, pos);
  				data.push_back(atoi(ss.c_str()));
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
  				data.push_back(atoi(ss.c_str()));
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
	    			data.push_back(atoi(ss.c_str()));
	    		}
				}
      }
      s.clear();
  	}
  }
  ifile.close();
}

void Read1D(string path, int index, int axis, vector<double>& data) {
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
  				data.push_back(atof(ss.c_str()));
  				s.erase(0, pos + 1);
  				ss.clear();
  				pos = s.find_first_of(',', 0);
  			}
  			pos = s.find_first_of('\n', 0);
  			if (pos != 0) {
  				ss = s.substr(0, pos);
  				data.push_back(atof(ss.c_str()));
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
  				data.push_back(atof(ss.c_str()));
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
	    			data.push_back(atof(ss.c_str()));
	    		}
				}
      }
      s.clear();
  	}
  }
  ifile.close();
}

void Print2D(string path, string mode, vector<vector<int> >& data) {
	ofstream ofile;
	if (mode == "app") {
		ofile.open(path.c_str(), ios::app);
	} else if (mode == "trunc") {
		ofile.open(path.c_str());
	}
	for (vector<vector<int> >::iterator it = data.begin(); it != data.end(); it++) {
		for (vector<int>::iterator itt = it->begin(); itt != it->end(); itt++) {
			ofile << *itt << ',';
		}
		ofile << '\n';
	}
	ofile.close();
}

void Print2D(string path, string mode, vector<vector<double> >& data) {
	ofstream ofile;
	if (mode == "app") {
		ofile.open(path.c_str(), ios::app);
	} else if (mode == "trunc") {
		ofile.open(path.c_str());
	}
	for (vector<vector<double> >::iterator it = data.begin(); it != data.end(); it++) {
		for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt++) {
			ofile << setprecision(15) << *itt << ',';
		}
		ofile << '\n';
	}
	ofile.close();
}

void Print2D(string path, string mode, vector<vector<bool> >& data) {
	ofstream ofile;
	if (mode == "app") {
		ofile.open(path.c_str(), ios::app);
	} else if (mode == "trunc") {
		ofile.open(path.c_str());
	}
	for (vector<vector<bool> >::iterator it = data.begin(); it != data.end(); it++) {
		for (vector<bool>::iterator itt = it->begin(); itt != it->end(); itt++) {
			ofile << *itt << ',';
		}
		ofile << '\n';
	}
	ofile.close();
}

void Print1D(string path, string mode, int axis, vector<int>& data) {
	ofstream ofile;
	if (mode == "app") {
		ofile.open(path.c_str(), ios::app);
	} else if (mode == "trunc") {
		ofile.open(path.c_str());
	}
	if (axis == 0) {
		for (vector<int>::iterator it = data.begin(); it != data.end(); it++) {
			ofile << *it << ',';
		}
		ofile << '\n';
	} else {
		for (vector<int>::iterator it = data.begin(); it != data.end(); it++) {
			ofile << *it << '\n';
		}
	}
	ofile.close();
}

void Print1D(string path, string mode, int axis, vector<double>& data) {
	ofstream ofile;
	if (mode == "app") {
		ofile.open(path.c_str(), ios::app);
	} else if (mode == "trunc") {
		ofile.open(path.c_str());
	}
	if (axis == 0) {
		for (vector<double>::iterator it = data.begin(); it != data.end(); it++) {
			ofile << setprecision(15) << *it << ',';
		}
		ofile << '\n';
	} else {
		for (vector<double>::iterator it = data.begin(); it != data.end(); it++) {
			ofile << setprecision(15) << *it << '\n';
		}
	}
	ofile.close();
}
