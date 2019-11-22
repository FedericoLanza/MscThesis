#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()
#include <random>
#include <vector>
#include <algorithm>

using namespace std;

unsigned int i;

void printvector(short int T, vector<int> v, ofstream &file ){
    for(unsigned int i=0; i<T; i++){
        file << v[i] << "\n";
    }
}


struct Xmin {
	double value;
	int xp;
};



int main(int argc, char *argv[])
{
	short int T = pow(2,10);  // size of the system; should be even
	short int dx = 2; //initial distance of the polymers; should be even
	unsigned int trials = 25*pow(10,4); // number of couple of polymers to be generated for the estimation of the non-crossing probability

	double N = 100000; // low energy polymers won't grow where we set energy to N

	//initializing the treshold array
	int ntresh = 0;
	for (i=0; i<T/2; i++){
		ntresh += i+1+dx/2;
	}
	ntresh = 2*ntresh + T/2+1+dx/2;
	vector<double> treshold(ntresh,0.0); //treshold array

	int k, time, x
	char avoid;
	vector<char> vavoid(trials,0.0); //the vector will count the number of pairs of polymer that avoid each other

	vector<vector<signed char>> lattice1(T+1, vector<signed char>(T+3,0)); //steps lattice for the left polymer
	vector<vector<signed char>> lattice2(T+1, vector<signed char>(T+3,0)); //steps lattice for the right polymer

	vector<double> E1(T+3,0.0); //energy of the g.s.
	vector<double> E2(T+3,0.0);
	Xmin Emin;
	vector<double> values(2,0.0);

	vector<int> Polym1(T+1,0); //polymer coordinates of the left g.s. (1)
	vector<int> Polym2(T+1,0); //polymer coordinates of the right g.s. (2)
	int itoto1, itoto2;
	signed char idelta1, idelta2;

	random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    //uniform_real_distribution<> udis(0.0, 1.0);	
	normal_distribution<double> ndis(0.0,1.0);

	for (unsigned int seed=0; seed<trials; seed++){

		//filling the treshold array with random values
		for (i=0; i<ntresh; i++){
			treshold[i]=ndis(gen);
		}
		k=0;

		//Dijkstra algorithm
		for (i=0; i<E1.size(); i++){
			E1[i] = N; E2[i] = N; // low energy polymers won't grow where we set energy to N=10^6
		}
		E1[start] = 0; E2[start] = 0;
    		for (time=1; time<=T/2; time++){
			for (x=T/2+1-time; x<=T/2+1+time; x=x+2){
				//left
				values = {E1[x-1], E1[x+1]};
				Emin.value = *min_element(values.begin(),values.end()); //take the minimum among the n.n.
				Emin.xp = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				if (Emin.xp == 0){
					Emin.xp = -1;
				}
        			else if (Emin.xp == 1){
        				Emin.xp = 1;
				}
                		lattice1[time][x] = Emin.xp;
				E1[x] = Emin.value+treshold[k];
				//right
				values = {E2[x-1], E2[x+1]};
				Emin.value = *min_element(values.begin(),values.end()); //take the minimum among the n.n.
				Emin.xp = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				if (Emin.xp == 0){
					Emin.xp = -1;
				}
        			else if (Emin.xp == 1){
        				Emin.xp = 1;
				}
                		lattice2[time][x] = Emin.xp;
				E2[x] = Emin.value+treshold[k+dx/2];
				k++;
			}
			k = k + dx/2;
		}
		for (time=T/2+1; time<=T; time++){
			for (x=time-T/2+1; x<=-time+3*T/2+1; x=x+2){
				//left
				values = {E1[x-1], E1[x+1]};
				Emin.value = *min_element(values.begin(),values.end()); //take the minimum among the n.n.
				Emin.xp = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				if (Emin.xp == 0){
					Emin.xp = -1;
				}
        			else if (Emin.xp == 1){
        				Emin.xp = 1;
				}
                		lattice1[time][x] = Emin.xp;
				E1[x] = Emin.value+treshold[k];
				//right
				values = {E2[x-1], E2[x+1]};
				Emin.value = *min_element(values.begin(),values.end()); //take the minimum among the n.n.
				Emin.xp = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				if (Emin.xp == 0){
					Emin.xp = -1;
				}
        			else if (Emin.xp == 1){
        				Emin.xp = 1;
				}
                		lattice2[time][x] = Emin.xp;
				E2[x] = Emin.value+treshold[k+dx/2];
				k++;
			}
			k = k + dx/2;
		}
		//cout << "k " << k << endl;
		//building of the g.s. polymers
    	 	itoto1 = T/2+1; //fix the arrival 
		itoto2 = T/2+1;
    		Polym1[T] = itoto1;
		Polym2[T] = itoto2 + dx;
		idelta1 = 0;
		idelta2 = 0;
		avoid = 0;
		time=1;
		//left
        	idelta1 = lattice1[T-time+1][itoto1]; //does the left polymer goes left or right?
        	itoto1 += idelta1;
        	Polym1[T-time] = itoto1; //position of the left polymer at time T-1-t
		//right
		idelta2 = lattice2[T-time+1][itoto2]; //does the right polymer goes left or right?
        	itoto2 += idelta2;
        	Polym2[T-time] = itoto2 + dx; //position of the right polymer at time T-1-t
    		while (Polym1[T-time]!=Polym2[T-time] && time!=T){
			time++;
			//left
        		idelta1 = lattice1[T-time+1][itoto1]; //does the left polymer goes left or right?
        		itoto1 += idelta1;
        		Polym1[T-time] = itoto1; //position of the left polymer at time T-1-t
			//right
			idelta2 = lattice2[T-time+1][itoto2]; //does the right polymer goes left or right?
        		itoto2 += idelta2;
        		Polym2[T-time] = itoto2 + dx; //position of the right polymer at time T-1-t
		}
		if (time==T){
			avoid = 1;
		}
		vavoid[seed] = avoid;
	}

	ofstream fileres ("pavoidnoise2D_results.txt");
    printvector(T, vavoid, fileres);	

	return 0;
}

