//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2017-03-09 16:55:42
//	Description: test program for get-config.h and get-config.cpp
//***************
#include "../get-config.h"  
#include <fstream>
#include <iostream>  

using namespace std;

int main() {
    string m_sPath = "./config.ini";  
    map<string,string> m_mapConfig;  
    ReadConfig(m_sPath,m_mapConfig);  
    PrintConfig(m_mapConfig);
    return 0;  
}  