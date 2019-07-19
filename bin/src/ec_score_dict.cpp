#ifndef EC_SCORE_DICT_CPP
#define EC_SCORE_DICT_CPP

// This is code to read ec score file, and store all the data in an ECSDict
// ECSDict : a dictionary-like data structure


#include <fstream>
#include <sstream>
//#include <iostream>
#include <string>
#include <unordered_map>
using namespace std;

// ECSDict[strandid][registration] = score
typedef unordered_map<int, unordered_map<int, double> > ECSDict;

void get_ec_dict(string fn,ECSDict& dict){
	string line;
	ifstream fin(fn);
	if(fin.fail()){
		cerr << "Error in opening file " << fn << endl;
		exit(1);
	}
	int strandi = 0;
	while(getline(fin,line)){
		dict[strandi] = unordered_map<int, double>();
		stringstream ss(line);
		int regnum;
		ss >> regnum;
		int reg;
		double score;
		for(int i=0; i<regnum; i++){
			ss >> reg >> score;
			dict[strandi][reg]=score;
		}
		strandi++;
	}
}

#endif//EC_SCORE_DICT_CPP
