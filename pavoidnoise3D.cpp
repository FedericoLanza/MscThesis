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


void printvector(short int T, vector<int> v, ofstream &file ){
    for(unsigned int i=0; i<T; i++){
        file << v[i] << "\n";
    }
}


int main(int argc, char *argv[])
{
	short int T = pow(2,10);  // size of the system; should be even
	short int dx = 4; //initial distance of the polymers; should be even
	short int a = 2; //distance under which the contact between the two polymer occur; should be even and > dx
	unsigned int trials = 25*pow(10,4); // number of couple of polymers to be generated for the estimation of the non-crossing probability

	double N = 100000; // low energy polymers won't grow where we set energy to N

	int time, y, n, x, i, j, k, shift
	char avoid;
	vector<char> vavoid(trials,0.0); //the vector will count the number of pairs of polymer that avoid each other
	double norma;
	vector<double> values(4,0.0); //energies of the n.n.
	double Emin; //minimum energy
	int EminIndex; //its position in 'values' vector
	signed char xmin, ymin;

	//steps lattice for the left polymer
	vector<vector<vector<signed char>>> xlattice1(T+1, vector<vector<signed char>>(T+3, vector<signed char>(T+3, 0))); 
	vector<vector<vector<signed char>>> ylattice1(T+1, vector<vector<signed char>>(T+3, vector<signed char>(T+3, 0)));
 	//steps lattice for the right polymer
	vector<vector<vector<signed char>>> xlattice2(T+1, vector<vector<signed char>>(T+3, vector<signed char>(T+3, 0)));
	vector<vector<vector<signed char>>> ylattice2(T+1, vector<vector<signed char>>(T+3, vector<signed char>(T+3, 0)));

	vector<vector<double>> E1(T+3, vector<double>(T+3)); //energy of the left polymer
	vector<vector<double>> E2(T+3, vector<double>(T+3)); //energy of the right polymer

	vector<short int> xPolym1(T+1,0); //coordinates of the left polymer
	vector<short int> yPolym1(T+1,0);
	vector<short int> xPolym2(T+1,0); //coordinates of the right polymer
	vector<short int> yPolym2(T+1,0);
	short int xitoto1, yitoto1, xitoto2, yitoto2;
	unsigned int seed, avoid;
	signed char xidelta1, yidelta1, xidelta2, yidelta2;
		
	//initializing the treshold array
	int ntresh = 0;
	for (i=1; i<T/2; i++){
		ntresh += pow(i+1+dx/2,2);
	}
	int triangle = 0;
	for (i=1; i<=dx/2; i++){
		triangle += i;
	}
	ntresh = 2*ntresh+pow(T/2+1+dx/2,2)+pow(dx/2,2)-2*T*triangle+2+dx;
	//cout << "ntresh: " << ntresh << endl;
	vector<double> treshold(ntresh,0.0); //treshold array

	random_device rd;  //Will be used to obtain a seed for the random number engine
    	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    	//uniform_real_distribution<> udis(0.0, 1.0);
	normal_distribution<double> ndis(0.0, 1.0);
	
	//ofstream filepart;

	for (seed=0; seed<trials; seed++){

		//restoring the energy
		for (i=0; i<T+3; i++){
			for (j=0; j<T+3; j++){
				E1[i][j] = N;
				E2[i][j] = N;
			}
		}
    		E1[T/2+1][T/2+1] = 0;
    		E2[T/2+1][T/2+1] = 0;

		//filling the treshold array with random values
		for (i=0; i<ntresh; i++){
			treshold[i]=ndis(gen);
		}
		k=0;

		//Dijkstra algorithm
    		for (time=1; time<=T/2; time++){
			for (y=T/2+1-time; y<=T/2+1+time; y++){
				for (x = T/2+1-time+abs(T/2+1-y); x <= T/2+1+time-abs(T/2+1-y); x=x+2){					
					//left
					values = {E1[y-1][x], E1[y+1][x], E1[y][x-1], E1[y][x+1]};
                			EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the min position
					Emin = *min_element(values.begin(), values.end()); //take the min among the n.n.
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
                			xlattice1[time][y][x] = xmin;
					ylattice1[time][y][x] = ymin;
                			E1[y][x] = Emin+treshold[k];
					//right
					values = {E2[y-1][x], E2[y+1][x], E2[y][x-1], E2[y][x+1]};
                			EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the min position
					Emin = *min_element(values.begin(), values.end()); //take the min among the n.n.
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
                			xlattice2[time][y][x] = xmin;
					ylattice2[time][y][x] = ymin;
                			E2[y][x] = Emin+treshold[k+dx/2];

					k++;
				}
				k=k+dx/2;
			}
		}
		for (time=T/2+1; time<=T; time++){
			for (y=time-T/2+1; y<=-time+3*T/2+1; y++){
				for (x = time-T/2+1+abs(T/2+1-y); x <= -time+3*T/2+1-abs(T/2+1-y); x=x+2){	
					//left
					values = {E1[y-1][x], E1[y+1][x], E1[y][x-1], E1[y][x+1]};
                			EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the min position
					Emin = *min_element(values.begin(), values.end()); //take the min among the n.n.
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
                			xlattice1[time][y][x] = xmin;
					ylattice1[time][y][x] = ymin;
                			E1[y][x] = Emin+treshold[k];
					//right
					values = {E2[y-1][x], E2[y+1][x], E2[y][x-1], E2[y][x+1]};
                			EminIndex = min_element(values.begin(),values.end()) - values.begin(); //take the min position
					Emin = *min_element(values.begin(), values.end()); //take the min among the n.n.
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
                			xlattice2[time][y][x] = xmin;
					ylattice2[time][y][x] = ymin;
                			E2[y][x] = Emin+treshold[k+dx/2];

					k++;
				}
				k=k+dx/2;
			}
		}
		//cout << "k: " << k << endl;
		//building of the g.s. polymer
    	 	xitoto1 = T/2+1;
		yitoto1 = T/2+1;
    	 	xitoto2 = T/2+1;
		yitoto2 = T/2+1;
    		xPolym1[T] = xitoto1; //fix the arrival at the same coordinate as of the beginning (i.e. (x,y)=0)
		yPolym1[T] = yitoto1;
    		xPolym2[T] = xitoto2 + dx; //fix the arrival at the same coordinate as of the beginning (i.e. (x,y)=0)
		yPolym2[T] = yitoto2;
		xidelta1 = 0;
		yidelta1 = 0;
		xidelta2 = 0;
		yidelta2 = 0;
		avoid = 0;
		time = 1;
		//left
        	xidelta1 = xlattice1[T-time+1][yitoto1][xitoto1]; //does the polymer goes left, right, up or down?
		yidelta1 = ylattice1[T-time+1][yitoto1][xitoto1];
        	xitoto1 += xidelta1;
		yitoto1 += yidelta1;
        	xPolym1[T-time] = xitoto1; //position of the left polymer at time T-t
		yPolym1[T-time] = yitoto1;
		//right
		xidelta2 = xlattice2[T-time+1][yitoto2][xitoto2]; //does the polymer goes left, right, up or down?
		yidelta2 = ylattice2[T-time+1][yitoto2][xitoto2];
        	xitoto2 += xidelta2;
		yitoto2 += yidelta2;
        	xPolym2[T-time] = xitoto2 + dx; //position of the right polymer at time T-t
		yPolym2[T-time] = yitoto2;
		norma = (xPolym1[T-time]-xPolym2[T-time])*(xPolym1[T-time]-xPolym2[T-time]) + (yPolym1[T-time]-yPolym2[T-time])*(yPolym1[T-time]-yPolym2[T-time]);
    		while (sqrt( norma ) > a && time != T){
			time++;
			//left
        		xidelta1 = xlattice1[T-time+1][yitoto1][xitoto1]; //does the polymer goes left, right, up or down?
			yidelta1 = ylattice1[T-time+1][yitoto1][xitoto1];
        		xitoto1 += xidelta1;
			yitoto1 += yidelta1;
        		xPolym1[T-time] = xitoto1; //position of the left polymer at time T-t
			yPolym1[T-time] = yitoto1;
			//right
			xidelta2 = xlattice2[T-time+1][yitoto2][xitoto2]; //does the polymer goes left, right, up or down?
			yidelta2 = ylattice2[T-time+1][yitoto2][xitoto2];
        		xitoto2 += xidelta2;
			yitoto2 += yidelta2;
        		xPolym2[T-time] = xitoto2 + dx; //position of the right polymer at time T-t
			yPolym2[T-time] = yitoto2;
			norma = (xPolym1[T-time]-xPolym2[T-time])*(xPolym1[T-time]-xPolym2[T-time]) + (yPolym1[T-time]-yPolym2[T-time])*(yPolym1[T-time]-yPolym2[T-time]);
		}
		if (time==T){
			avoid = 1;
		}
		vavoid[seed] = avoid;
	}

	ofstream fileres ("pavoidnoise3D_results.txt");
    printvector(T, vavoid, fileres);	

	return 0;
}




