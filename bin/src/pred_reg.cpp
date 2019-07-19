#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include "ec_score_dict.cpp"
#include "anyoption.h"
#include "anyoption.cpp"

using namespace std;

enum StrandOrientation {PERI2EXTRA,EXTRA2PERI};
enum HBondPattern {SWV,VWS};

// load odds data
void load_odds(const string fn, vector<double>& arr, double default_val, bool needlog=true);
// calculate pairing energy of a hairpin
double pairing(HBondPattern hbp, int len, StrandOrientation o1, vector<int> &strand1, vector<int> &strand2) ;
// get pairid for a given aa pair
int aa_pair(int a, int b);
// determine hbond pattern of two neighboring strands
HBondPattern patterning(int extra, int orientation1);
int getres(int pdbi, int seqid);


const static int AABoundary = 19;
const static int InvalidAA = 20;
const static double nref = 8.5; // avg loop length of ori 25 prot


vector<double> strongs(210,0), vdws(210,0), weaks(210,0), CoreIn(21,0), CoreOut(21,0), PeriIn(21,0), PeriOut(21,0), ExtraIn(21,0), ExtraOut(21,0);
vector<int> res;

int main(int argc,char*argv[]) {
	AnyOption *opt = new AnyOption();
	opt->setOption("class"    );
	opt->setOption("res"      );
	opt->setOption("strand"   );
	opt->setOption("seqcov"   );
	opt->setOption("strong"   );
	opt->setOption("weak"     );
	opt->setOption("vdw"      );
	opt->setOption("extraout" );
	opt->setOption("extrain"  );
	opt->setOption("coreout"  );
	opt->setOption("corein"   );
	opt->setOption("periout"  );
	opt->setOption("periin"   );

	opt->processCommandArgs(argc, argv);

	int cls = 0;
	string resfn, strandfn, seqcovfn, strongfn, weakfn, vdwfn, extraoutfn, extrainfn, coreoutfn, coreinfn, perioutfn, periinfn;

	if(opt->getValue("class"    )) { cls = atoi(opt->getValue("class")); }

	if(opt->getValue("res"      )) { resfn      = opt->getValue("res"      ); }
	if(opt->getValue("seqcov"   )) { seqcovfn   = opt->getValue("seqcov"   ); }
	if(opt->getValue("strand"   )) { strandfn   = opt->getValue("strand"   ); }

	if(opt->getValue("strong"   )) { strongfn   = opt->getValue("strong"   ); } else{ strongfn   = "odds/all.strong.odds"; }
	if(opt->getValue("weak"     )) { weakfn     = opt->getValue("weak"     ); } else{ weakfn     = "odds/all.weak.odds";   }
	if(opt->getValue("vdw"      )) { vdwfn      = opt->getValue("vdw"      ); } else{ vdwfn      = "odds/all.vdw.odds";    }
	if(opt->getValue("extraout" )) { extraoutfn = opt->getValue("extraout" ); } else{ extraoutfn = "odds/ExtraOut.odds";   }
	if(opt->getValue("extrain"  )) { extrainfn  = opt->getValue("extrain"  ); } else{ extrainfn  = "odds/ExtraIn.odds";    }
	if(opt->getValue("coreout"  )) { coreoutfn  = opt->getValue("coreout"  ); } else{ coreoutfn  = "odds/CoreOut.odds";    }
	if(opt->getValue("corein"   )) { coreinfn   = opt->getValue("corein"   ); } else{ coreinfn   = "odds/CoreIn.odds";     }
	if(opt->getValue("periout"  )) { perioutfn  = opt->getValue("periout"  ); } else{ perioutfn  = "odds/PeriOut.odds";    }
	if(opt->getValue("periin"   )) { periinfn   = opt->getValue("periin"   ); } else{ periinfn   = "odds/PeriIn.odds";     }


	double wec, penub, penneg, wstrong, wvdw, wweak; // weights
	switch(cls){
		case 0:
			wec=1; penub=0.245; penneg=0.05; wstrong=0.026; wvdw=0.038; wweak=0.036;
			break;
		case 1:
			wec=1; penub=0.45; penneg=0.12; wstrong=0.055; wvdw=0.1; wweak=0.075;
			break;
		case 2:
			wec=1; penub=0.052; penneg=0.074; wstrong=0; wvdw=0.082; wweak=0.006;
			break;
		case 3:
			wec=1; penub=0.29; penneg=0.1; wstrong=0.045; wvdw=0.02; wweak=0.024;
			break;
		case 4:
			wec=1; penub=0.11; penneg=0.135; wstrong=0.045; wvdw=0.024; wweak=0.014;
			break;
		default:
			//TODO error
			break;
	}


	// load res data
	res.push_back(InvalidAA); // index 0
	ifstream fin(resfn);
	int aares;
	while(fin >> aares){
		aares = aares>AABoundary?InvalidAA:aares;
		res.push_back(aares);
	}
	fin.close();


	// read ec score
	ECSDict dict;
	get_ec_dict(seqcovfn,dict);


	vector<int> peris, extras;
	// read strands
	fin.open(strandfn);
	int end1, end2;
	int strand_num = 0;
	while(fin >> end1 >> end2){
		if(strand_num%2==0){
			peris.push_back(end1);
			extras.push_back(end2);
		}
		else{
			peris.push_back(end2);
			extras.push_back(end1);
		}
		strand_num++;
	}


	// load pair odds
	load_odds(strongfn, strongs, -1.72);
	load_odds(vdwfn, vdws, -1.72);
	load_odds(weakfn, weaks, -1.72);
	for(int i=0; i<210; i++) {
		strongs[i] *= wstrong;
		vdws[i] *= wvdw;
		weaks[i] *= wweak;
	}

	// load single body odds
	load_odds( extraoutfn, ExtraOut, -3.9 );
	load_odds( extrainfn,  ExtraIn, -3.9 );
	load_odds( coreoutfn, CoreOut, -3.9 );
	load_odds( coreinfn,  CoreIn, -3.9 );
	load_odds( perioutfn, PeriOut, -3.9 );
	load_odds( periinfn,  PeriIn, -3.9 );


	// loop for each strand for current pdb
	for(int strandi = 0; strandi < strand_num; strandi++) {
		int strandj = (strandi+1) % strand_num;
		vector<int> newstrand1, newstrand2;
		StrandOrientation orientation1;
		if(strandi%2==0){
			orientation1 = PERI2EXTRA;
			// make strands
			newstrand1 = vector<int> (res.begin()+peris[strandi], res.begin()+extras[strandi]+1);
			newstrand2 = vector<int> (res.begin()+extras[strandj], res.begin()+peris[strandj]+1);
			reverse(newstrand1.begin(),newstrand1.end());
		}
		else{
			orientation1 = EXTRA2PERI;
			// make strands
			newstrand1 = vector<int> (res.begin()+extras[strandi], res.begin()+peris[strandi]+1);
			newstrand2 = vector<int> (res.begin()+peris[strandj], res.begin()+extras[strandj]+1);
			reverse(newstrand2.begin(),newstrand2.end());
		}

		// determine pattern
		HBondPattern hbond_pattern = patterning(extras[strandi], orientation1);

		// make tmp strand2
		vector<int> newstrand2tmp(newstrand1.size());

		double maxscore = -1000;
		int predreg  =0;
		// enumerate registration
		for(int regoffset=-10; regoffset <= 6; regoffset++) { 
			for(int i=0; i < newstrand1.size(); i++) {
				if(i+regoffset < 0 || i+regoffset >= newstrand2.size())
					newstrand2tmp[i] = InvalidAA;
				else
					newstrand2tmp[i] = newstrand2[i+regoffset];
			}

			int currreg = newstrand2.size()-newstrand1.size()-regoffset;
			double ecscore = dict[strandi][currreg];
			double negative_reg = currreg>=0 ? 0 : currreg;
			double currscore = pairing(hbond_pattern, newstrand1.size(), orientation1, newstrand1, newstrand2tmp) // hbond
								- penub*log((abs(regoffset) + nref)/nref) // penalty for unbonded res
								+ penneg * negative_reg // panalty for neg reg
								+ wec * ecscore; // ec score

			cout << setprecision(3) << currreg << ":" << currscore << " ";

			if(currscore > maxscore) {
				maxscore = currscore;
				predreg = currreg;
			}
		}
		//cout << " " << strandi << " " << predreg << endl;
		cout << endl;
		//cout << "### " << predreg << endl;
	}// strand loop
	return 0;
}

/// 
int aa_pair(int aa1, int aa2) {
	if(aa1 > AABoundary || aa2 > AABoundary ) { return -1; }
	if(aa1 <= aa2) { return (aa1 * 20 + aa2 - (aa1*(aa1+1)/2)); }
	else { return (aa2 * 20 + aa1 - (aa2*(aa2+1)/2)); }
}

double pairing(HBondPattern hbp, int len, StrandOrientation o1, vector<int> &strand1, vector<int> &strand2) {
	double sum = 0.0;
	vector<int> &tmpstrand1 = (o1==PERI2EXTRA) ? strand2 : strand1;
	vector<int> &tmpstrand2 = (o1==PERI2EXTRA) ? strand1 : strand2;
	if(hbp == SWV) {
		for(int i = 0 ; i < len-1 ; i++) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i+1]) != -1) {
				sum += weaks[aa_pair(tmpstrand1[i], tmpstrand2[i+1])];
			}
		}
		for(int i = 0 ; i < len ; i = i+2) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i]) != -1) {
				sum += strongs[aa_pair(tmpstrand1[i], tmpstrand2[i])];
			}
		}
		for(int i = 1 ; i < len ; i = i+2) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i]) != -1) {
				// will: gly no vdw? diff from pairing()
				if( tmpstrand1[i] != 7 || tmpstrand2[i] != 7)
					sum += vdws[aa_pair(tmpstrand1[i], tmpstrand2[i])]; 
			}
		}
	}
	else {
		for(int i = 0 ; i < len-1 ; i++) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i+1]) != -1) {
				sum += weaks[aa_pair(tmpstrand1[i], tmpstrand2[i+1])];
			}
		}
		for(int i = 0 ; i < len ; i = i+2) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i]) != -1) {
				if( tmpstrand1[i] != 7 || tmpstrand2[i] != 7)
					sum += vdws[aa_pair(tmpstrand1[i], tmpstrand2[i])]; 
			}
		}
		for(int i = 1 ; i < len ; i = i+2) {
			if(aa_pair(tmpstrand1[i], tmpstrand2[i]) != -1) {
				sum += strongs[aa_pair(tmpstrand1[i], tmpstrand2[i])];
			}
		}	
	}
	return sum;
}



void load_odds(const string fn, vector<double>& arr, double default_val, bool needlog){
	ifstream fin(fn.c_str());
	if(fin.fail()) {
		cerr << "error  opening file " << fn <<endl;
		exit(1);
	}
	int i = 0;
	while(fin >> arr[i]){
		if(arr[i]==0){ arr[i]=default_val; }
		else{ arr[i] = needlog ? (double)log(arr[i]) : (double)arr[i]; }
		i++;
	}
	fin.close();
}


int getres(int seqid){
	if( seqid<1 || seqid>res.size()-1 ){
		return InvalidAA;
	}
	return res[seqid];
}


HBondPattern patterning(int extra, int orientation1){
	double SumEven, SumOdd;
	if(orientation1 == EXTRA2PERI) {
		SumEven = ExtraIn[getres(extra)]+ExtraOut[getres(extra+1)]+CoreIn[getres(extra+2)]+CoreOut[getres(extra+3)]+CoreIn[getres(extra+4)]+CoreOut[getres(extra+5)]+CoreIn[getres(extra+6)]+PeriOut[getres(extra+7)]+PeriIn[getres(extra+8)];
		SumOdd = ExtraOut[getres(extra)]+ExtraIn[getres(extra+1)]+CoreOut[getres(extra+2)]+CoreIn[getres(extra+3)]+CoreOut[getres(extra+4)]+CoreIn[getres(extra+5)]+CoreOut[getres(extra+6)]+PeriIn[getres(extra+7)]+PeriOut[getres(extra+8)];
	}
	else {
		SumOdd = ExtraIn[getres(extra)]+ExtraOut[getres(extra-1)]+CoreIn[getres(extra-2)]+CoreOut[getres(extra-3)]+CoreIn[getres(extra-4)]+CoreOut[getres(extra-5)]+CoreIn[getres(extra-6)]+PeriOut[getres(extra-7)]+PeriIn[getres(extra-8)];
		SumEven = ExtraOut[getres(extra)]+ExtraIn[getres(extra-1)]+CoreOut[getres(extra-2)]+CoreIn[getres(extra-3)]+CoreOut[getres(extra-4)]+CoreIn[getres(extra-5)]+CoreOut[getres(extra-6)]+PeriIn[getres(extra-7)]+PeriOut[getres(extra-8)];
	}
	if(SumOdd > SumEven){
		return SWV;
	}
	else{
		return VWS;
	}
}

