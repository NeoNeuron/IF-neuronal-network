#include<string>
#include<vector>
using namespace std;

//	Read 1D data from *.txt files; Data type: (double);
//	Return: none;
void ReadData(string filename, vector<double> & data);

// Return mean value of elements of x;
double Mean(vector<double>& x);

// Return standard deviation of elements of x;
double Std(vector<double>& x);

// Pearson's product-moment coefficient: Linear correlated coefficient;
double LC(vector<int>& raster, vector<double>& lfp);

void TDLC(vector<int>& raster, vector<double>& lfp, int negative_time_delay, int positive_time_delay, vector<double>& tdlc);
