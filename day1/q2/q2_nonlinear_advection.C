#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

int main()
{
	int i, itr;
	int scheme;
	double dx, dt, x;
	double lambda;

	vector<double> f,fnew;

	dx = 0.001;
	dt = 0.00001;
	double maxT = 1.0;
	int n = 1000;
	double itr_max = maxT / dt;

	lambda = dt / dx;

	for(i = 0;i<n+1;i++){
		f.push_back(0);
		fnew.push_back(0);
	}

	for(i = 0;i<n;i++){
		if(i>=100 && i<=300){ 
			//below commented is for sine wave
			// if(i>=100 && i<=(2*M_PI/dx+100)){      
			// x = (i-100) * dx;
			//f[i] = (10 * x - 4) * (10 * x - 4) * ( 6 - 10 * x) * (6- 10 * x);
			f[i] = 1.0;
			//     cout<<i<<" "<<f[i]<<endl;
		}
		else
			f[i] = 0.0;

		}

		ofstream myFn;
		myFn.open("heaviside.dat");
		for(i = 0;i<n;i++)myFn<<i*dx<<" "<<f[i]<<endl;


		cout<<"Choose your favorite scheme"<<endl;
		cout<<"1=FTCS"<<" "<<"2=FTBS"<<" "<<"3=Lax"<<" "<<"4=Lax-Wendroff"<<endl;
		cout<<"5=Beam-warming"<<" "<<"6=Mac-Cormak"<<endl;
		cin>>scheme;


		if(scheme == 1){//FTCS
			for (itr=1; itr<=itr_max; itr++){
				
				// Write FTCS Scheme here
				
											}
						}

		if(scheme == 2){//FTBS
			for (itr=1; itr<=itr_max; itr++){
				
				// Write FTBS Scheme here
				
											}
						}


		if(scheme == 3){//Lax-Friedrichs
			for (itr=1; itr<=itr_max; itr++){
				
				// Write Lax-Friedrichs Scheme here
				
											}
						}

		if(scheme == 4){// Lax - wendroff
			for (itr=1; itr<=itr_max; itr++){
				
				// Write Lax - wendroff Scheme here
				
											}
						}

		if(scheme == 5){// Beam-Warming
			for (itr=1; itr<=itr_max; itr++){
				
				// Write Beam-Warming Scheme here
				
											}
						}

		if(scheme == 6){//Mac-cormak 
			vector<double> fp;
			for(i = 0;i<n+1;i++) fp.push_back(0);
			for (itr=1; itr<=itr_max; itr++){
				
				// Write Mac-cormak Scheme here
				
											}
						}


		std::ostringstream fileName;
		fileName << "scheme_NLC_" << scheme << ".dat";
		std::ofstream myFile(fileName.str().c_str(),std::ios_base::binary);
		if (myFile.is_open()){
			for (i=0; i<n; i++)
				myFile<<(i*dx)<<" "<<f[i]<<std::endl;
		}
		myFile.close();

		return 0;

	}

