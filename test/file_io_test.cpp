#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

int main() {
  // Create test sequence;
  bool data[10];
  for (int i = 0; i < 10; i ++) {
    data[i] = false;
    cout << data[i] << endl;
  }
  // output data into bin file;
  ofstream bfile;
  ofstream tfile;
  bfile.open("test.bin", ios::binary);
  tfile.open("test.csv");
  if (bfile) {
    // for (vector<vector<double> >::iterator it = data.begin(); it != data.end(); it ++) {
      // for (vector<double>::iterator itt = it->begin(); itt != it->end(); itt ++) {
      //   bfile.write((char*)&*itt, sizeof(*itt));
      // }
      // bfile.write((char*)&*it, it->size()*8);
    // }
    bfile.write((char*)&data, 10*sizeof(bool));
    // cout << sizeof(data) << endl;
    // cout << bfile.tellp() << endl;
  } else cout << "Fail to open target file." << endl;
  bfile.close();
  // if (tfile) {
  //   for (vector<vector<double> >::iterator it = data.begin(); it != data.end(); it ++) {
  //     tfile << setprecision(15) << *(it->begin()) << ',' << *(it->begin() + 1) << '\n';
  //   }
  // } else cout << "Fail to open target file." << endl;
  // tfile.close();

  ifstream ibfile;
  ibfile.open("test.bin", ios::binary);
  // vector<double> buffer(2);
  bool buffer[10];
  // ibfile.seekg(32, ibfile.beg);
  ibfile.read((char*)&buffer, 10*sizeof(bool));
  for (int i = 0; i < 10; i ++) {
    cout << buffer[i] << ',';
  }
  cout << endl;
  cout << sizeof(buffer) << endl;
  cout << sizeof(signed char) << endl;
  cout << sizeof(short int) << endl;
  
  ibfile.close();
  return 0;
}
