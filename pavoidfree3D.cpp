#include <iostream>     // std::cout
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <algorithm>	// std::random_shuffle
#include <complex>	// std::complex

using namespace std;

void printvector(short int T, vector<int> v, ofstream &file ){
    for(unsigned int i=0; i<T; i++){
        file << v[i] << "\n";
    }
}

int round_tailcut(double x, double T){
	int n;
	if (x<0){
		n=0;
	}
	else if(x>T/2){
		n=T/2;
	}
	else{
		n=round(x);
	}
	return n;
}

int main () {
	short int T = pow(2,10);  // size of the system; should be even
	short int dx = 4; //initial distance of the polymers; should be even
	short int a = 2; //distance under which the contact between the two polymer occur; should be even and > dx
	unsigned int trials = 25*pow(10,4); // number of couple of polymers to be generated for the estimation of the non-crossing probability

	vector<char> navoid; navoid = 0.0;
	char avoid;
	int n, m, k, j, t;
	srand ( unsigned ( time(0) ) );
	vector<complex<double>> stepsa(T,(0,0));
	vector<complex<double>> stepsb(T,(0,0));
	vector<complex<double>> patha(T+1,(0,0));
	vector<complex<double>> pathb(T+1,(0,0));
	complex<double> uno(1,0);
	complex<double> muno(-1,0);
	complex<double> i(0,1);
	complex<double> mi(0,-1);
	patha[0] = (dx/2)*uno;
	pathb[0] = (dx/2)*muno;
	random_device rd;  // Will be used to obtain a seed for the random number engine
    	mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    	normal_distribution<double> ndis(T/4,(T/4+1)/sqrt(T+6));
	for (m=0; m<trials; m++){
		// generating the arrays of the steps
		n=round_tailcut(ndis(gen),T);
		for (j=0; j<2*n; j=j+2){
			stepsa[j] = uno;
			stepsa[j+1] = muno;
		}
		for (j=2*n; j<T; j=j+2){
			stepsa[j] = i;
			stepsa[j+1] = mi;
		}
		random_shuffle ( stepsa.begin(), stepsa.end() );
		n=round_tailcut(ndis(gen),T);
		for (j=0; j<2*n; j=j+2){
			stepsb[j] = uno;
			stepsb[j+1] = muno;
		}
		for (j=2*n; j<T; j=j+2){
			stepsb[j] = i;
			stepsb[j+1] = mi;
		}
		random_shuffle ( stepsb.begin(), stepsb.end() );
		// building the paths and verifying if they avoid
		t = 1;
		avoid = 0;
		patha[t] = patha[t-1] + stepsa[t-1];
		pathb[t] = pathb[t-1] + stepsb[t-1];
		while (sqrt(norm(patha[t] - pathb[t])) > a && t != T){
			t++;
			patha[t] = patha[t-1] + stepsa[t-1];
			pathb[t] = pathb[t-1] + stepsb[t-1];
			}
		if (t == T){
			avoid=1;
		}
		vavoid[m]=avoid;
	}

	ofstream fileres ("pavoidnoise3D_results.txt");
    printvector(T, vavoid, fileres);

	return 0;
}


