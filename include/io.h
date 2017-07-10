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
using namespace std;

//	Read 2 dimensional information;
	//	Read 2-D data from *.csv files; data type can be int or double;
	//	Return: none;
void Read2D(string path, vector<vector<int> >& data);
void Read2D(string path, vector<vector<double> >& data);

//	Read 1 dimensional (line/column) from *.csv file;
  //	STRING path: file path of *.csv file;
  //	INT index: the index of chosen line/column, from 0 to maximum - 1;
  //  INT axis: 0 for line, 1 for column;
  //	VECTOR<INT/DOUBLE> data: storage of output data;
  //	Return: none;
void Read1D(string path, int index, int axis, vector<int>& data);
void Read1D(string path, int index, int axis, vector<double>& data);

//  Print 2 dimensional data;
  //  Print 2-D data to *.csv files; data type can be int, double or bool;
  //  STRING mode: openmode for aiming file; "app" for append, "trunc" for any contents that existed in the file before it is open are discarded.
  //  Return: none;
void Print2D(string path, string mode, vector<vector<int> >& data);
void Print2D(string path, string mode, vector<vector<double> >& data);
void Print2D(string path, string mode, vector<vector<bool> >& data);

//	Print 1 dimensional (line/column) to *.csv file;
  //	STRING path: file path of *.csv file;
  //  STRING mode: openmode for aiming file; "app" for append, "trunc" for any contents that existed in the file before it is open are discarded.
  //  INT axis: 0 for line, 1 for column;
  //	VECTOR<INT/DOUBLE> data: storage of outputing data;
  //	Return: none;
void Print1D(string path, string mode, int axis, vector<int>& data);
void Print1D(string path, string mode, int axis, vector<double>& data);
#endif // _IO_H_
