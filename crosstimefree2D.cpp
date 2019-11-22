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
	int T = pow(2,8);  // size of the system; should be even
	int dx = 4; //initial distance of the polymers; should be even
	double N = 100000; // low energy polymers won't grow where we set energy to N
	unsigned int trials = 100000; // number of reseaches of g.s. + first excited polymers
	vector<double> timecross(trials,0.0); // array of energy gaps
	
	//initializing the treshold array
	int ntresh = 0;
	for (int i=1; i<T+1; i++){
		ntresh=ntresh+i+2;
	}
	vector<double> treshold(ntresh,0.0); //treshold array

	int i, k, time, x, xmin, ymin, EminIndex;
	double crosses, pavoid, navoid;
	navoid = 0.0;
	vector<double> values(2,0.0); //energies of the n.n.
	double Emin; //minimum energy

	vector<vector<signed char>> lattice1(T+1, vector<signed char>(2*T+1,0)); //steps lattice for the left polymer
	vector<vector<signed char>> lattice2(T+1, vector<signed char>(2*T+1,0)); //steps lattice for the right polymer

	vector<double> E1(2*T+1,0.0); //energy of the g.s.
	vector<double> E2(2*T+1,0.0);

	vector<int> Polym1(T+1,0); //polymer coordinates of the left g.s. (1)
	vector<int> Polym2(T+1,0); //polymer coordinates of the right g.s. (2)
	int itoto1, itoto2;
	signed char idelta1, idelta2;

	random_device rd;  //Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    //uniform_real_distribution<> udis(0.0, 1.0);	
	normal_distribution<double> ndis(0.0,1.0);

	for (unsigned int seed=0; seed<trials; seed++){
		
		//restoring the energy
		for (i=0; i<2*T+1; i++){
			E1[i] = N;
			E2[i] = N;
		}
    		E1[T] = 0;
		E2[T] = 0;
		


		for (i=0; i<ntresh; i++){
			treshold[i]=dis(gen);
		}
		k=0;

		//Dijkstra algorithm
    		for (time=1; time<T; time++){
			for (x=T-time; x<=T+time; x=x+2){
				//left
				values = {E1[x-1], E1[x+1]};
                		EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				Emin = minimo(values); //take the minimum among the n.n.
                		if (EminIndex == 0){
		            		xmin = -1;
				}
                		else if (EminIndex == 1){
                    			xmin = 1;
				}
                		lattice1[time][x] = xmin;
				E1[x] = Emin+treshold[k];
				//right
				values = {E2[x-1], E2[x+1]};
                		EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
				Emin = minimo(values); //take the minimum among the n.n.
                		if (EminIndex == 0){
		            		xmin = -1;
				}
                		else if (EminIndex == 1){
                    			xmin = 1;
				}
                		lattice2[time][x] = xmin;
				E2[x] = Emin+treshold[k + dx/2];
				k++;
			}
			k = k + dx/2;
		}
		time = T;
		for (x=T-time; x<=T+time; x=x+2){
			//left
			if (x==0){
				values = {N, E1[x+1]};
			}
			else if (x==2*T){
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
               		lattice1[time][x] = xmin;
               		E1[x] = Emin+treshold[k];
			//right
			if (x==0){
				values = {N, E2[x+1]};
			}
			else if (x==2*T){
				values = {E2[x-1], N};
			}
			else{
				values = {E2[x-1], E2[x+1]};
			}
               		EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the position of the min
			Emin = minimo(values); //take the minimum among the n.n.
               		if (EminIndex == 0){
			    	xmin = -1;
			}
               		else if (EminIndex == 1){
                	   	xmin = 1;
			}
               		lattice2[time][x] = xmin;
               		E2[x] = Emin+treshold[k + dx/2];

			k++;
		}

		//building of the g.s. polymers
    	 	itoto1 = T + dx/2;
		itoto2 = T - dx/2;
    		Polym1[T] = itoto1; //fix the arrival 
		Polym2[T] = itoto2 + dx; 
		idelta1 = 0;
		idelta2 = 0;
		crosses = 0;
    		time=1;
		//left
        	idelta1 = lattice1[T-time+1][itoto1]; //does the left polymer goes left or right?
        	itoto1 += idelta1;
        	Polym1[T-time] = itoto1; //position of the left polymer at time T-1-t
		//right
		idelta2 = lattice2[T-time+1][itoto2]; //does the right polymer goes left or right?
        	itoto2 += idelta2;
        	Polym2[T-time] = itoto2 + dx; //position of the right polymer at time T-1-t
    		while (Polym1[T-time]==Polym2[T-time] && time!=T){
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
		timecross[seed]=time-1;
	}
	cout << "T = " << T << ", dx = " << dx << ", trials = " << trials << endl;

	ofstream file ("polymdis2D.txt");
    	printPolym(Polym1, Polym2, T+1, file);
	file.close();

	ofstream file2 ("crosstime2D_T256.txt");
	printvec(timecross, file2);
	file2.close();

	return 0;
}
