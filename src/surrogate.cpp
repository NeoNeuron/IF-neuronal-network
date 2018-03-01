#include <iostream>
#include <algorithm>
#include <ctime>
#include "../include/surrogate.h"
using namespace std;

void Shuffle(vector<double>& data) {
  srand(time(NULL));
  random_shuffle(data.begin(), data.end());
}

void Shuffle2D1D(vector<vector<double> >& data, int axis){
  if (axis == 0) {
    time_t rseed = time(NULL);
    for (vector<vector<double> >::iterator it = data.begin(); it != data.end(); it ++) {
      srand(rseed);
      random_shuffle(it->begin(), it->end());
    }
  } else if (axis == 1) {
    srand(time(NULL));
    random_shuffle(data.begin(), data.end());
  }
}
