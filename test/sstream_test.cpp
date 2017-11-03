//  This program aims to test the function of stringstream class. It differentiate the usage of str("") and clear();
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
using namespace std;

int main() {
  ifstream ifile("tfile.txt");
  string s;
  istringstream ss;
  string buffer;
  vector<bool> boollist;
  int col_counter = 0;
  while (getline(ifile, s)) {
    cout << s << ' ';
    col_counter = 0;
    ss.clear();
    ss.str("");
    ss.str(s);
    // istringstream ss(s);
    cout << ss.str() << ' ';
    getline(ss,buffer,',');
    cout << buffer << ' ' << ss.str() << '\n';
    // while (getline(ss,buffer,',')) {
      // boollist.push_back(atoi(buffer.c_str()));
      // cout << buffer << endl;
    // }

  }
  ifile.close();
  return 0;
}
