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
#include <string>
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
  istringstream s_input;
  while (getline(ifile, s)) {
    add_T.clear();
    s_input.clear();
    s_input.str("");
    s_input.str(s);
    while (getline(s_input, buffer, ',')) {
      if (buffer.find('.') == string::npos) {
        add_T.push_back(atoi(buffer.c_str()));
      } else add_T.push_back(atof(buffer.c_str()));
    }
    data.push_back(add_T);
  }
  ifile.close();
}

template <class T> void Read2DBin(string path, vector<vector<T> >& data) {
  // open data file;
  ifstream ifile(path.c_str(), ios::binary);
  // Read the shape of 2D data array:
  size_t shape[2];
  ifile.read((char*)&shape, 8);
  // Clear the inputing vector;
  data.clear();
  data.resize(shape[0], vector<T>(shape[1]));
  T buffer;
  size_t size = sizeof(buffer);
  for (size_t i = 0; i < shape[0]; i ++) {
    for (size_t j = 0; j < shape[1]; j ++) {
      ifile.read((char*)&buffer, size);
      data[i][j] = buffer;
    }
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
    istringstream s_input;
    while (getline(ifile, s)) {
      if (getline_counter == index) {
        s_input.clear();
        s_input.str("");
        s_input.str(s);
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
    istringstream s_input;
    while (getline(ifile, s)) {
      column_counter = 0;
      s_input.clear();
      s_input.str("");
      s_input.str(s);
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

template <class T> void Read1DBin(string path, vector<T>& data, int index, int axis) {
  // open data file;
  ifstream ifile(path.c_str(), ios::binary);
  // Read the shape of 2D data array:
  size_t shape[2];
  ifile.read((char*)&shape, 8);
  // clear the data vector;
  data.clear();
  T buffer;
  size_t size = sizeof(T);
  if (axis == 0) {
    size_t prefix = index * shape[1];
    ifile.seekg(prefix + 8, ifile.beg);
    for (int i = 0; i < shape[1]; i ++) {
      ifile.read((char*)&buffer, size);
      data.push_back(buffer);
    }
  } else {
    ifile.seekg(index + 8, ifile.beg);
    for (size_t i = 0; i < shape[0]; i ++) {
      ifile.read((char*)&buffer, size);
      data.push_back(buffer);
      ifile.seekg(shape[1] - 1, ifile.cur);
    }
  }
  ifile.close();
}

template <class T> void Print2D(string path, vector<vector<T> >& data, string mode) {
  ofstream ofile;
  if (mode == "app") ofile.open(path.c_str(), ios::app);
  else if (mode == "trunc") ofile.open(path.c_str());
  ostringstream s_out;
  for (typename vector<vector<T> >::iterator it = data.begin(); it != data.end(); it++) {
    s_out.str("");
    for (typename vector<T>::iterator itt = it->begin(); itt != it->end(); itt++) {
      s_out << setprecision(15) << *itt << ',';
    }
    s_out << '\n';
    ofile << s_out.str();
  }
  ofile.close();
}

template <class T> void Print2DBin(string path, vector<vector<T> >& data, string mode) {
  ofstream ofile;
  if (mode == "app") ofile.open(path.c_str(), ios::binary|ios::app);
  else if (mode == "trunc") {
    ofile.open(path.c_str(), ios::binary);
    size_t buffer;
    buffer = data.size();
    ofile.write((char*)&buffer, 4);
    buffer = data.begin() -> size();
    ofile.write((char*)&buffer, 4);
  }
  size_t size = sizeof(T);
  T temp;
  for (typename vector<vector<T> >::iterator it = data.begin(); it != data.end(); it++) {
    for (typename vector<T>::iterator itt = it -> begin(); itt != it -> end(); itt ++) {
      temp = *itt;
      ofile.write((char*)&temp, size);
    }
  }
  ofile.close();
}

template <class T> void Print1D(string path, vector<T>& data, string mode, int axis) {
  ofstream ofile;
  if (mode == "app") ofile.open(path.c_str(), ios::app);
  else if (mode == "trunc") ofile.open(path.c_str());
  ostringstream s_out;
  if (axis == 0) {
    for (typename vector<T>::iterator it = data.begin(); it != data.end(); it++) {
      s_out << setprecision(15) << *it << ',';
    }
    s_out << '\n';
    ofile << s_out.str();
  } else {
    for (typename vector<T>::iterator it = data.begin(); it != data.end(); it++) {
      s_out << setprecision(15) << *it << '\n';
    }
    ofile << s_out.str();
  }
  ofile.close();
}

template <class T> void Print1DBin(string path, vector<T>& data, string mode) {
  ofstream ofile;
  if (mode == "app") ofile.open(path.c_str(), ios::binary|ios::app);
  else if (mode == "trunc") {
    ofile.open(path.c_str(), ios::binary);
    size_t buffer;
    buffer = 1;
    ofile.write((char*)&buffer, 4);
    buffer = data.size();
    ofile.write((char*)&buffer, 4);
  }
  size_t size = sizeof(T);
  T temp;
  for (typename vector<T>::iterator it = data.begin(); it != data.end(); it ++) {
    temp = *it;
    ofile.write((char*)&temp, size);
  }
  ofile.close();
}
