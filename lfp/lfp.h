//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-08 21:32:46
//	Description: Define model of local field potential(LFP);
//***************
#ifndef _DATA_ANALYSIS_LFP_H_
#define _DATA_ANALYSIS_LFP_H_

#include <string>
#include <vector>

using namespace std;

//	Read specific line from *.txt file;
//	STRING filename: file name of txt file;
//	INT line_index: the index of chosen lines, from 0 to maximum - 1;
//	VECTOR<INT> output: storage of output data;
//	Return: none;
void ReadLine(string filename, int line_index, vector<int> &output);

//	Read specific line from *.txt file;
//	STRING filename: file name of txt file;
//	INT line_index: the index of chosen lines, from 0 to maximum - 1;
//	VECTOR<DOUBLE> output: storage of output data;
//	Return: none;
void ReadLine(string filename, int line_index, vector<double> &output);

//	Read specific Column from *.txt file;
//	STRING filename: file name of txt file;
//	INT column_index: the index of chosen lines, from 0 to maximum - 1;
//	INT num_colunm: total number of column contained in target file; nature number;
//	VECTOR<INT> output: storage of output data;
//	Return: none;
void ReadColumn(string filename, int column_index, int num_column, vector<int> &output);

//	Read specific lines from *.txt file;
//	STRING filename: file name of txt file;
//	VECTOR<VECTOR<DOUBLE> > data: storage for the lines of data;
//	VECTOR<INT> line_index: List of lines that choosen;
//	Return: none;
void ReadLines(string filename, vector<int> &line_index, vector<vector<double> > &data);

//	ChooseNeurons:
//	Choose the specific portion of neurons in post-network;
//	


//	Local field potential model [version 0.10]
//	Description: point current source model without sptial distribution;
//	DOUBLE* t_range: time period used in calculation, with unit ms, include the last point while not the first point
//	STRING potential_file: membrane potential;
//	STRING excitatory_conductance_file;
//	STRING inhibitory_conductance_file;
//	VECTOR<DOUBLE> lfp: local field potential data;
//	Return: none;
void LFP(double* t_range, vector<int> & neuron_list, string potential_filename, string excitatory_conductance_filename, string inhibitory_conductance_filename, vector<double> &lfp);

// 	Output LFP;
//	Description: output LFP data to a given file;
//	STRING filename: output file name;
//	Return: none;
void OutputLFP(vector<double> &lfp, string filename);

// 	Output spike train;
//	Description: output the spike train of chosen neuron within chosen time range; spiking time points are rearranged which starts from 0;
//	DOUBLE* t_range: time period used in calculation, with unit ms;
//	STRING filename: output file name;
//	Return: none;
void OutputSpikeTrain(double* t_range, vector<double> &spikes, string filename);

#endif // _DATA_ANALYSIS_LFP_H_