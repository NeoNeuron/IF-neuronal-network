//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-08-29
//	Description: Library for Data I/O functions in project;
//***************
#ifndef _IO_H_
#define _IO_H_

#include <string>
#include <vector>
using namespace std;

//	Read 2 dimensional information;
	//	Read 2-D data from *.csv files; data type can be int or double;
  //	VECTOR<VECTOR<T> > data: container of data;
	//	STRING path: path of target file;
	//	Return: none;
template <class T> void Read2D(string path, vector<vector<T> >& data);
template <class T> void Read2DBin(string path, vector<vector<T> >& data);
//	Read 1 dimensional (line/column) from *.csv file;
  //	STRING path: file path of *.csv file;
  //	VECTOR<INT/DOUBLE> data: storage of output data;
  //	INT index: the index of chosen line/column, from 0 to maximum - 1;
  //  INT axis: 0 for line, 1 for column;
  //	Return: none;
template <class T> void Read1D(string path, vector<T>& data, int index, int axis);
template <class T> void Read1DBin(string path, vector<T>& data, int index, int axis);

//  Print 2 dimensional data;
  //  Print 2-D data to *.csv files; data type can be int, double or bool;
  //	VECTOR<VECTOR<T> > data: container of data;
  //  STRING mode: openmode for aiming file; "app" for append, "trunc" for any contents that existed in the file before it is open are discarded.
  //  Return: none;
template <class T> void Print2D(string path, vector<vector<T> >& data, string mode);
template <class T> void Print2DBin(string path, vector<vector<T> >& data, string mode);

//	Print 1 dimensional (line/column) to *.csv file;
  //	STRING path: file path of *.csv file;
  //	VECTOR<INT/DOUBLE> data: storage of outputing data;
  //  STRING mode: openmode for aiming file; "app" for append, "trunc" for any contents that existed in the file before it is open are discarded.
  //  INT axis: 0 for line, 1 for column;
  //	Return: none;
template <class T> void Print1D(string path, vector<T>& data, string mode, int axis);
template <class T> void Print1DBin(string path, vector<T>& data, string mode);

#include "io.hpp"

#endif // _IO_H_
