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

void printPolym( vector<int> v1, vector<int> v2, int T, ofstream &file ){
    for(unsigned int i=0; i<T; i++){
        file << v1[i] << " " << v2[i] << " " << i << "\n";
    }
}

void printPolymD( vector<int> v1, vector<int> v2, vector<int> v3, vector<int> v4, int T, ofstream &file ){
    for(unsigned int i=0; i<T; i++){
        file << v1[i] << " " << v2[i] << " " << v3[i] << " " << v4[i] << " " << i << "\n";
    }
}

void printvec( vector<double> v, ofstream &file ){
    for(unsigned int i=0; i<v.size(); i++){
        file << v[i] << "\n";
    }
}

double minimo(vector<double> v){
	double min; min=v[0];
	for(int i=0; i<v.size(); i++){
		if (v[i]<min){
			min=v[i];
		}
	}
	return min;
}


int main(int argc, char *argv[])
{
	int T = 32;  // size of the system; should be odd (and >=3)
	double N = 100000; // low energy polymers won't grow where we set energy to N
	unsigned int trials = 1000000; // number of reseaches of g.s. + first excited polymers
	vector<double> dE(trials,0.0); // array of energy gaps
	
	//initializing the treshold array
	int ntresh = 0;
	for (int i=1; i<T+1; i++){
		ntresh=ntresh+i+2;
	}
	vector<double> treshold(ntresh,0.0); //treshold array

	int time, x, i, k;

	vector<double> values(2,0.0); //energies of the n.n.
	double Emin; //minimum energy
	int EminIndex; //its position in 'values' vector
	int xmin, ymin;

	vector<vector<signed char>> lattice(T+1, vector<signed char>(T+1)); //steps lattice

	vector<double> E(T+1,0.0); //energy of the g.s.
	vector<double> E1(T+1,0.0); //energy of the first excited

	vector<int> Polym(T+1,T/2); //polymer coordinates of the g.s.
	vector<int> Polym1(T+1,T/2); //polymer coordinates of the first excited
	int itoto;
	signed char idelta;

	//filling the treshold array with random values
	random_device rd;  //Will be used to obtain a seed for the random number engine
    	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    	uniform_real_distribution<> dis(0.0, 1.0);

	for (unsigned int seed=0; seed<trials; seed++){

		for (i=0; i<ntresh; i++){
			treshold[i]=dis(gen);
		}

		//G.S. POLYMER
		k=0;
		//restoring the energy
		for (i=0; i<T+1; i++){
			E[i] = N;
		}
		E[T/2-1] = treshold[k];
		k++;
		E[T/2+1] = treshold[k];
		k++;
		E[T/2+1] = N;

		//Dijkstra algorithm
    		for (time=2; time<T/2; time++){
			for (x=T/2-time; x<=T/2+time; x=x+2){
				values = {E[x-1], E[x+1]};
                		EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				Emin = minimo(values); //take the minimum among the n.n.
                		if (EminIndex == 0){
		            		xmin = -1;
				}
                		else if (EminIndex == 1){
                    			xmin = 1;
				}
                		lattice[time][x] = xmin;
                		E[x] = Emin+treshold[k];
				k++;
			}
		}
		time=T/2;
		for (x=time-T/2; x<=-time+3*T/2; x=x+2){
			if (x == 0){
			values = {N, E[x+1]};
			}
			else if (x == T){
				values = {E[x-1], N};
			}
			else{
				values = {E[x-1], E[x+1]};
			}
               		EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
			Emin = minimo(values); //take the minimum among the n.n.
               		if (EminIndex == 0){
			    	xmin = -1;
			}
               		else if (EminIndex == 1){
                	   	xmin = 1;
			}
               		lattice[time][x] = xmin;
               		E[x] = Emin+treshold[k];
			k++;
		}
		for (time=T/2+1; time<T; time++){
			for (x=time-T/2; x<=-time+3*T/2; x=x+2){
				values = {E[x-1], E[x+1]};
               			EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				Emin = minimo(values); //take the minimum among the n.n.
               			if (EminIndex == 0){
				    	xmin = -1;
				}
               			else if (EminIndex == 1){
                		   	xmin = 1;
				}
               			lattice[time][x] = xmin;
               			E[x] = Emin+treshold[k];
				k++;
			}
		}
		//building of the g.s. polymer
    	 	itoto = T/2-1;
    		Polym[T-1] = itoto; //fix the arrival at the same coordinate as of the beginning
		idelta = 0;
    		for (time=2; time<T; time++){
        		idelta = lattice[T+1-time][itoto]; //does the polymer goes left, right, up or down?
        		itoto+=idelta;
        		Polym[T-time] = itoto; //position of the polymer at time T-1-t
		}


		//Ist EXCITED POLYMER
		//restore energy
		for (i=0; i<T+1; i++){
			E1[i] = N;
		}
		k=0;
		E1[T/2-1] = treshold[k];
		k++;
		E1[T/2+1] = treshold[k];
		k++;
		E1[T/2-1] = N;

		//Dijkstra algorithm
    		for (time=2; time<T/2; time++){
			for (x=T/2-time; x<=T/2+time; x=x+2){
				values = {E1[x-1], E1[x+1]};
                		EminIndex = min_element(values.begin(),values.end()) - values.begin();
				Emin = minimo(values);
        	        	if (EminIndex == 0){
			        	xmin = -1;
				}
        	        	else if (EminIndex == 1){
        	            		xmin = 1;
				}
        	        	lattice[time][x] = xmin;
        	        	E1[x] = Emin+treshold[k];
				k++;
			}
			E1[Polym[time]] = N;
		}
		time=T/2;
		for (x=time-T/2; x<=-time+3*T/2; x=x+2){
			if (x == 0){
				values = {N, E1[x+1]};
			}
			else if (x == T){
				values = {E1[x-1], N};
			}
			else{
				values = {E1[x-1], E1[x+1]};
			}
               		EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
			Emin = minimo(values); //take the minimum among the n.n.
               		if (EminIndex == 0){
			    	xmin = -1;
			}
               		else if (EminIndex == 1){
                	   	xmin = 1;
			}
               		lattice[time][x] = xmin;
               		E1[x] = Emin+treshold[k];
			k++;
			E1[Polym[time]] = N;
		}
		for (time=T/2+1; time<T; time++){
			for (x=time-T/2; x<=-time+3*T/2; x=x+2){
				values = {E1[x-1], E1[x+1]};
               			EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				Emin = minimo(values); //take the minimum among the n.n.
               			if (EminIndex == 0){
				    	xmin = -1;
				}
               			else if (EminIndex == 1){
                		   	xmin = 1;
				}
               			lattice[time][x] = xmin;
               			E1[x] = Emin+treshold[k];
				k++;
			}
			E1[Polym[time]] = N;
		}
    		//building of the first excited polymer
    	 	itoto = T/2+1;
    		Polym1[T-1] = itoto; //fix the arrival at the same coordinate as of the beginning (i.e. x=0)
		idelta = 0;
    		for (time=2; time<T; time++){
        		idelta = lattice[T+1-time][itoto]; //does the polymer goes left or right?
        		itoto+=idelta;
        		Polym1[T-time] = itoto; //position of the polymer at time T-1-t
		}
		
		dE[seed] = E1[T/2+1]-E[T/2-1];

	}

	ofstream file ("dE_T32.txt");
	printvec(dE, file);
	return 0;
}

