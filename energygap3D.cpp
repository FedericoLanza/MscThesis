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

void printdE( vector<double> v, ofstream &file ){
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
	int T = 65;  // size of the system; should be odd (and >=3)
	double N = 100000; // low energy polymers won't grow where we set energy to N
	unsigned int trials = 1000000; // number of reseaches of g.s. + first excited polymers

	vector<double> dE(trials,0.0); // array of energy gaps
	//initializing the random treshold array
	vector<int> floor(T-1,0);
	floor[0] = 4;
	for (unsigned int i=1; i<floor.size(); i++){
		if (i % 2 == 0){
			floor[i] = floor[i-1]+(i/2)*4+3;
		}
		else{
			floor[i] = floor[i-1]+(i+1)*2+1;
		}
	}
	
	int ntresh = accumulate(floor.begin(), floor.end(), 0); //sum of all floors
	vector<double> treshold(ntresh,0.0); //treshold array

	int time, y, n, x, i, j, k;

	vector<double> values(4,0.0); //energies of the n.n.
	double Emin; //minimum energy
	int EminIndex; //its position in 'values' vector
	int xmin, ymin;
	signed char xlattice[T][2*T-1][2*T-1]; //steps lattice
	signed char ylattice[T][2*T-1][2*T-1];
	for (time=0; time<T; time++){		
		for (i=0; i<2*T-1; i++){
			for (j=0; j<2*T-1; j++){
				xlattice[time][i][j] = 0;
				ylattice[time][i][j] = 0;
			}
		}
	}
	double E[2*T-1][2*T-1]; //energy of the g.s.
	double E1[2*T-1][2*T-1]; //energy of the first excited

	vector<int> xPolym(T,0); //polymer coordinates of the g.s.
	vector<int> yPolym(T,0);
	vector<int> xPolym1(T,0); //polymer coordinates of the first excited
	vector<int> yPolym1(T,0);
	int xitoto;
	int yitoto;
	signed char xidelta;
	signed char yidelta;

	for (unsigned int seed=0; seed<trials; seed++){

		
		//FIRST POLYMER

		//restoring the energy
		for (i=0; i<2*T-1; i++){
			for (j=0; j<2*T-1; j++){
				E[i][j] = N;
			}
		}
    		E[T-1][T-1] = 0;

		//filling the random treshold array
		random_device rd;  //Will be used to obtain a seed for the random number engine
    		mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    		uniform_real_distribution<> dis(0.0, 1.0);
		for (i=0; i<ntresh; i++){
			treshold[i]=dis(gen);
		}
		k=0;

		//Dijkstra algorithm
    		for (time=1; time<=T-2; time++){
			for (y=T-1-time; y<=T-1+time; y++){
				for (n=0; n<=time-abs(T-1-y); n++){
					x = T-1+time-abs(T-1-y)-2*n;
					values = {E[y-1][x], E[y+1][x], E[y][x-1], E[y][x+1]};
                			EminIndex = min_element(values.begin(),values.end()) - values.begin();
					Emin = minimo(values);
                			if (EminIndex == 0){
						xmin = 0;
		            			ymin = -1;
					}
                			else if (EminIndex == 1){
                    				xmin = 0;
                    				ymin = 1;
					}
					else if (EminIndex == 2){
						xmin = -1;
		            			ymin = 0;
					}
                			else if (EminIndex == 3){
                    				xmin = 1;
                    				ymin = 0;
					}
                			xlattice[time][y][x] = xmin;
					ylattice[time][y][x] = ymin;
                			E[y][x] = Emin+treshold[k];
					k++;
				}
			}
		}
		time = T-1;
		for (y=T-1-time; y<=T-1+time; y++){
			for (n=0; n<=time-abs(T-1-y); n++){
				x = T-1+time-abs(T-1-y)-2*n;
				if (y==0 && x==T-1){
					values = {N, E[y+1][x], E[y][x-1], E[y][x+1]};
				}
				else if (y==T-1 && x==0){
					values = {E[y-1][x], E[y+1][x], N, E[y][x+1]};
				}
				else if (y==T-1 && x==2*T-2){
					values = {E[y-1][x], E[y+1][x], E[y][x-1], N};
				}
				else if (y==2*T-2 && x==T-1){
					values = {E[y-1][x], N, E[y][x-1], E[y][x+1]};
				}
				else{
					values = {E[y-1][x], E[y+1][x], E[y][x-1], E[y][x+1]};
				}
               			EminIndex = min_element(values.begin(),values.end()) - values.begin();
				Emin = minimo(values);
               			if (EminIndex == 0){
					xmin = 0;
			           	ymin = -1;
				}
               			else if (EminIndex == 1){
                		   	xmin = 0;
                		   	ymin = 1;
				}
				else if (EminIndex == 2){
					xmin = -1;
			           	ymin = 0;
				}
               			else if (EminIndex == 3){
                		   	xmin = 1;
                		   	ymin = 0;
				}
               			xlattice[time][y][x] = xmin;
				ylattice[time][y][x] = ymin;
               			E[y][x] = Emin+treshold[k];
				k++;
			}
		}

		//building of the g.s. polymer
    	 	xitoto = T-1;
		yitoto = T-1;
    		xPolym[T-1] = xitoto; //fix the arrival at the same coordinate as of the beginning (i.e. (x,y)=0)
		yPolym[T-1] = yitoto;
		xidelta = 0;
		yidelta = 0;
    		for (time=1; time<=T-1; time++){
        		xidelta = xlattice[T-time][yitoto][xitoto]; //does the polymer goes left, right, up or down?
			yidelta = ylattice[T-time][yitoto][xitoto];
        		xitoto+=xidelta;
			yitoto+=yidelta;
        		xPolym[T-1-time] = xitoto; //position of the polymer at time T-1-t
			yPolym[T-1-time] = yitoto;
		}


		//SECOND POLYMER
		//restore some arrays
		for (i=0; i<=2*T-2; i++){
			for (j=0; j<=2*T-2; j++){
				E1[i][j] = N;
			}
		}
    		E1[T-1][T-1] = 0;
		k=0;

		//Dijkstra algorithm
    		for (time=1; time<=T-2; time++){
			for (y=T-1-time; y<=T-1+time; y++){
				for (n=0; n<=time-abs(T-1-y); n++){
					x = T-1+time-abs(T-1-y)-2*n;
					values = {E1[y-1][x], E1[y+1][x], E1[y][x-1], E1[y][x+1]};
                			EminIndex = min_element(values.begin(),values.end()) - values.begin();
					Emin = minimo(values);
                			if (EminIndex == 0){
						xmin = 0;
		            			ymin = -1;
					}
                			else if (EminIndex == 1){
                    				xmin = 0;
                    				ymin = 1;
					}
					else if (EminIndex == 2){
						xmin = -1;
		            			ymin = 0;
					}
                			else if (EminIndex == 3){
                    				xmin = 1;
                    				ymin = 0;
					}
                			xlattice[time][y][x] = xmin;
					ylattice[time][y][x] = ymin;
                			E1[y][x] = Emin+treshold[k];
					E1[yPolym[time]][xPolym[time]] = N;
					k++;
				}
			}
		}
		time = T-1;
		for (y=T-1-time; y<=T-1+time; y++){
			for (n=0; n<=time-abs(T-1-y); n++){
				x = T-1+time-abs(T-1-y)-2*n;
				if (y==0 && x==T-1){
					values = {N, E[y+1][x], E1[y][x-1], E1[y][x+1]};
				}
				else if (y==T-1 && x==0){
					values = {E1[y-1][x], E1[y+1][x], N, E1[y][x+1]};
				}
				else if (y==T-1 && x==2*T-2){
					values = {E1[y-1][x], E1[y+1][x], E1[y][x-1], N};
				}
				else if (y==2*T-2 && x==T-1){
					values = {E1[y-1][x], N, E1[y][x-1], E[y][x+1]};
				}
				else{
					values = {E1[y-1][x], E1[y+1][x], E1[y][x-1], E1[y][x+1]};
				}
               			EminIndex = min_element(values.begin(),values.end()) - values.begin();
				Emin = minimo(values);
               			if (EminIndex == 0){
					xmin = 0;
		           		ymin = -1;
				}
               			else if (EminIndex == 1){
                   			xmin = 0;
                   			ymin = 1;
				}
				else if (EminIndex == 2){
					xmin = -1;
		           		ymin = 0;
				}
               			else if (EminIndex == 3){
                   			xmin = 1;
                   			ymin = 0;
				}
               			xlattice[time][y][x] = xmin;
				ylattice[time][y][x] = ymin;
               			E1[y][x] = Emin+treshold[k];
				k++;
			}
		}
		
		dE[seed] = E1[T-1][T-1]-E[T-1][T-1];

	}

	ofstream file ("dE.txt");
	printdE(dE, file);
	return 0;
}




