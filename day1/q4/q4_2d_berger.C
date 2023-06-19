#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include<bits/stdc++.h>
using namespace std;
int main()
{
	int i,j, itr;
	//int scheme;
	double dx, dy,dt, x;
	double lambda;


	dx = 0.01;
	dy = dx;
	dt = 0.001;
	double maxT = 1.0;
	int n = 1000;
	double itr_max = maxT / dt;
	double nu = 0.01;

	lambda = dt / dx;
	static  double u[1001][1001],unew[1001][1001],v[1001][1001],vnew[1001][1001];
	static double up [1001][1001] ;
	static double vp [1001][1001] ;


	for(int i=0;i<1001;i++)
	{
		for(int j=0;j<1001;j++)
		{
			u[i][j]=0.0;
			unew[i][j]=0.0;
			v[i][j]=0.0;
			vnew[i][j]=0.0;
			up[i][j]=0.0;
			vp[i][j]=0.0;
		}
	}

	for(i = 100;i<401;i++)
	{
		for (j=0; j<n; j++)
		{
			//below commented is for sine wave
			// if(i>=100 && i<=(2*M_PI/dx+100)){      
			// x = (i-100) * dx;
			//f[i] = (10 * x - 4) * (10 * x - 4) * ( 6 - 10 * x) * (6- 10 * x);
			u[i][j] = 1.0;
			v[i][j]= 1.0;
			//     cout<<i<<" "<<f[i]<<endl;
		}

		/* else
		   u[i][j] = 0.0;
		   v[i][j] = 0.0;*/

		}
		ofstream myFn;
		myFn.open("heaviside_2d.dat");
		for(i = 0;i<n;i++)
			for (j=1; j<n; j++)
				myFn<<i*dx<<" "<<i*dy<<" "<<u[i][j]<<" "<<v[i][j]<<endl;


		//Mac-cormack Scheme

		for (itr=1; itr<=itr_max; itr++){
			for (i=1; i<n; i++)
			{
				for (j=1; j<n; j++)
				{
					//...........
				}
			}

			for (i=1; i<n; i++)
			{
				for (j=1; j<n; j++)
				{
					//.........
				}}

			for (i=1; i<n; i++)
			{
				for (j=1; j<n; j++)
				{
					//......
				}}
		}


		std::ostringstream fileName;
		fileName <<"burger_2d.dat";
		std::ofstream myFile(fileName.str().c_str(),std::ios_base::binary);
		if (myFile.is_open()){
			for (i=0; i<n; i++)
				for (j=0; j<n; j++)
					myFile<<(i*dx)<<" "<<(j*dy)<<" "<<u[i][j]<<" "<<v[i][j]<<" "<<sqrt(pow(u[i][j],2) +pow(v[i][j],2))<< std::endl;
		}
		myFile.close();

		return 0;


	}
