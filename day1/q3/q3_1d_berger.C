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

	dx = 0.01;
	dt = 0.000001;
	double maxT = 1.0;
	int n = 1000;
	double itr_max = maxT / dt;
	double nu = 0.01;

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

		/*      
			cout<<"Choose your favorite scheme"<<endl;
			cout<<"1=FTCS"<<" "<<"2=FTBS"<<" "<<"3=Lax"<<" "<<"4=Lax-Wendroff"<<endl;
			cout<<"5=Beam-warming"<<" "<<"6=Mac-Cormak"<<endl;
			cin>>scheme;


			if(scheme == 1){//FTCS
			for (itr=1; itr<=itr_max; itr++){
			for (i=1; i<n; i++)
			fnew[i] = f[i] - 0.5 * lambda*f[i] * (f[i+1]-f[i-1]);

			for (i=1; i<n; i++)f[i] = fnew[i];
			f[n] = f[1];     f[0] = f[n-1];}}

			if(scheme == 2){//FTBS
			for (itr=1; itr<=itr_max; itr++){
			for (i=2; i<=n; i++)
			fnew[i] = f[i] - lambda *f[i]* (f[i]-f[i-1]);

			for (i=2; i<=n; i++)f[i] = fnew[i];
			f[1] = f[n];f[0] = f[n-1];}}


			if(scheme == 3){//Lax-Friedrichs
			for (itr=1; itr<=itr_max; itr++){
			for (i=1; i<n; i++)
			fnew[i] = 0.5 * (f[i+1] + f[i-1] - lambda *f[i]* (f[i+1]-f[i-1]));

			for (i=1; i<n; i++)f[i] = fnew[i];
			f[n] = f[1];     f[0] = f[n-1];}}

			if(scheme == 4){// Lax - wendroff
			for (itr=1; itr<=itr_max; itr++){
			for (i=1; i<n; i++)
			fnew[i] = f[i] -0.5 * lambda *f[i]* (f[i+1]-f[i-1]) + 0.5 * lambda * lambda *f[i]*f[i]* (f[i+1]-2*f[i]+f[i-1]);

			for (i=1; i<n; i++)f[i] = fnew[i];
			f[n] = f[1];     f[0] = f[n-1];}}

			if(scheme == 5){// Beam-Warming
			for (itr=1; itr<=itr_max; itr++){
			for (i=2; i<=n; i++)
			fnew[i] = (1 - 1.5 * lambda*f[i]  + 0.5 * lambda * lambda*f[i]*f[i]) * f[i] + lambda *f[i]* (2 - lambda*f[i]) * f[i-1] + 0.5 * lambda*f[i] * (lambda*f[i] - 1) * f[i-2];

			for (i=2; i<=n; i++)f[i] = fnew[i];
			f[0] = f[n-2]; f[1] = f[n-1];}}
		 */

		//mac-cormack
		vector<double> fp;
		for(i = 0;i<n+1;i++) fp.push_back(0);
		for (itr=1; itr<=itr_max; itr++){
				
				// Write Mac-cormak Scheme here
				
										}


		std::ostringstream fileName;
		fileName << "burger_1d.dat";
		std::ofstream myFile(fileName.str().c_str(),std::ios_base::binary);
		if (myFile.is_open()){
			for (i=0; i<n; i++)
				myFile<<(i*dx)<<" "<<f[i]<<std::endl;
		}
		myFile.close();

		return 0;

	}

