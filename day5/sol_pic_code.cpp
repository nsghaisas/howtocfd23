#include <iostream>
#include <fstream>
#include<fftw3.h>
#include<string>
#include<cmath>
#define _USE_MATH_DEFINES
#include<random>
#include <iomanip>

using namespace std;
using std::setw;

int nrows, ncols;  //number of rows, number of columns
//Macros for real and img part
#define REAL 0
#define IMAG 1
int nppc, nc;      //number of particles, number of cells
double Lx;         //Length of the system
double dx;         //width of each cell
int v_0;           //initial velocity
double** ParArray;     // Array with particle position, charge and mass
double* velocity_x;     //velocity array
double* charge_density;   //charge density at grid points
double* E_field;
double* xgrid;
double* pcharge_density;
double* echarge_density;
int init();
int calc_charge_density();  //calculating charge density
int output();
int E_cal();
int deallocate();           // deallocating arrays
int half_velocity();
int particle_push();
double dt;                  // time-step  
double mi,mr;            // electron mass, proton mass, mass ratio
double ntime;                  // Total time
double* times;              // time array
int dout;                   // interval for output
//
double* rho_t;
int test_rho();
int t;

int main() {
   
    int ret = init(); //calling initialization function

    cout << "total time:  " << endl;
    cin >> ntime;
    cout << "Interval for output:  " <<endl;
    cin >> dout;
    dt = 0.001;
    int itime = ntime/dt;
    times = new double[itime];
    // Take a half-step to get v_(1/2)
    int rho = calc_charge_density();
    int ret_E = E_cal();
    int hf_v = half_velocity();
    //cause of routine to advance by half step-----to do
    for (t=0; t<itime;t++) {
         times[t]=t*dt;
         rho = calc_charge_density();  //calculating charge density
         ret_E = E_cal(); //calling Efield calculation
         int ptp=particle_push();//particle pusher
         if (t%dout==0) {
             int ret_o = output();
         }
	 cout<<times[t]<<endl;
    }


    int free = deallocate(); //deallocating arrays

    delete[] times;
    times = NULL;
   
    return 0;
}

int init() { //developer Sruti
    cout << "Initializing simulation" << endl;
    cout << "Number of particles per cell :" << endl;
    cin >> nppc;
    cout << "Number of cells: " << endl;
    cin >> nc;
    cout << "Length of the system: " << endl;
    cin >> Lx;
    cout << "Mass ratio: " << endl;
    cin >> mr;
    nrows = 3;
    ncols = 2*nppc*nc; //total number of particles
    dx = Lx/nc;
    v_0 = 0; 
    mi = 1;
    double dq=0.01;
    double kt=(2*M_PI)/Lx;

    default_random_engine rand;
    int a=0;
    int b=Lx;


    // Position array
    ParArray = new double*[nrows];
    for (int i=0; i<nrows; i++) {
         ParArray[i] = new double[ncols];
    }

    // Storing user input in the array
    for (int j = 0; j < ncols; ++j) {
         //if (j%2==0) {
         ParArray[0][j]=(j*Lx)/ncols + 0.000001;  // TODO position 
    }
    

    for (int j = 0; j < ncols; j++) {   //ions
         if (j%2==0) {
             ParArray[1][j]=1.0/nppc; //TODO charge [+1, (rho*dx)] //kdm changes introduced division by nppc
             ParArray[2][j]=1.0/nppc; //TODO mass [(1*dx)]
         }
         else {   //elctrons
             ParArray[1][j]=(-1.0 + dq*sin(kt*ParArray[0][j]))/nppc; //TODO charge [-1, (rho*dx)]
             ParArray[2][j]=(1.0/mr)/nppc; //TODO  mass [(dx/mass ratio)----electron]
         }
    }

    E_field = new double[nc+1];         // Electric field at grid points
    for (int j=0; j< nc+1; j++) {
         E_field[j]=0;
    }


    velocity_x = new double[ncols];      //velocity array
    for (int j=0; j<ncols; j++) {
	 velocity_x[j] =0;
    }

    xgrid = new double[nc+1];
    for (int j=0; j<nc+1; j++) {
         xgrid[j]=j*dx;
    }

    charge_density = new double[nc+1];     //charge density at grid points
    for (int j=0; j< nc+1; j++) {
	 charge_density[j]=0;
    }

    pcharge_density = new double[nc+1];     //proton charge density at grid points
    for (int j=0; j< nc+1; j++) {
         pcharge_density[j]=0;
    }

    echarge_density = new double[nc+1];     //electron charge density at grid points
    for (int j=0; j< nc+1; j++) {
         echarge_density[j]=0;
    }

    return 0;
}

int calc_charge_density() { // calculating charge density; charge density stored in array index 0 to nc

    //double rcharge[nc+1];     //charge density at grid points
    double pcharge[nc+1];     //proton charge density at grid points
    double echarge[nc+1];     //electron charge density at grid points


    for (int i=0; i<nc+1; i++) {
         //rcharge[i]=0;
	 pcharge[i]=0;
	 echarge[i]=0;
    }

    //TODO deposit charge of particles into charge density on the grid points
    
    for (int j = 0; j < ncols; j ++) {

         int leftgrid  = (ParArray[0][j]) / dx;
         int rightgrid = leftgrid + 1;
         
         double q_i  = ParArray[1][j] * ((xgrid[rightgrid] - ParArray[0][j]) / dx) / dx;//kdm divide by dx to get density
         double q_i1 = ParArray[1][j] * ((ParArray[0][j]  - xgrid[leftgrid]) / dx) / dx;
         
         if (j%2 == 0) { //ions
             pcharge[leftgrid]  += q_i;
             pcharge[rightgrid] += q_i1;
         }
         else {   //elctrons
             echarge[leftgrid]  += q_i;
             echarge[rightgrid] += q_i1;
         }    
    }

    double eghost_cell = echarge[nc]+echarge[0];
    double pghost_cell = pcharge[nc]+pcharge[0];
    echarge[0]=eghost_cell;
    echarge[nc]=eghost_cell;
    pcharge[0]=pghost_cell;
    pcharge[nc]=pghost_cell;

    for (int j=0;j<nc+1;j++) {
         charge_density[j]=pcharge[j]+echarge[j];
	 pcharge_density[j]=pcharge[j];
	 echarge_density[j]=echarge[j];
	 //cout << "ec:  "<<echarge_density[j] << endl;
    }

    return 0;
}

int deallocate() { //deallocating arrays
    for (int i=0; i<nrows; i++) {
        delete[] ParArray[i];
    }
    delete[] ParArray;
    ParArray = NULL;
    delete[] velocity_x;
    velocity_x=NULL;
    delete[] xgrid;
    xgrid=NULL;
    delete[] charge_density;
    charge_density=NULL;
    delete[] pcharge_density;
    pcharge_density=NULL;
    delete[] echarge_density;
    echarge_density=NULL;
    delete[] E_field;
    E_field=NULL;
    return 0;
}

int output() { //developer Johan
    string t_s= to_string(t);;
    cout << "Output" << endl;

    ofstream fout("output/output"+t_s+".txt");
    
    //Storing the output to a text file
    if(fout.is_open())
    {
        for (int i=0; i<nc+1; i++){
        fout<<std::scientific;
	if (i==0){
		fout<<"Grid position"<<"  "<<setw(20)<<"Proton charge density"<<"  "<<setw(20)<<"Electron charge density"<<"   "<<setw(17)<<"Electric Field"<<endl;
		fout<<xgrid[i]<<"  "<<setw(20)<<pcharge_density[i]<<"  "<<setw(20)<<echarge_density[i]<<"   "<<setw(20)<<E_field[i]<<endl;

	}
	else{
		fout<<xgrid[i]<<"  "<<setw(20)<<pcharge_density[i]<<"  "<<setw(20)<<echarge_density[i]<<"   "<<setw(20)<<E_field[i]<<endl;
	}}
    cout << "Success!" << endl;
    }
    else
    {
    cout << "File could not be opened." << endl;
    }
    return 0;
}

int E_cal() { //developer Johan
        double kx[nc+1];
	double delkx=(2*M_PI)/Lx;

        //input array
        fftw_complex rho_x[nc+1];     //Charge density at grid point in real space
        //Output array
        fftw_complex rho_k[nc+1];      //Fourier transform of charge density
        fftw_complex E_k[nc+1];    //Electric field in fourier space
	fftw_complex E_x[nc+1];    //Electric field in real space
	fftw_complex Ex_norm[nc+1]; //Normalised electric field in real space
        
	//fill the first array with charge density
        for (int i=0;i<nc+1;i++)
        {
                rho_x[i][REAL] = charge_density[i]; //charge_density[i];
                rho_x[i][IMAG] = 0;
        }
        //Plant the fft and execute it (of charge density)
        fftw_plan plan = fftw_plan_dft_1d(nc+1,rho_x,rho_k,FFTW_FORWARD,FFTW_ESTIMATE);
        fftw_execute(plan);
        //Do some cleaning
        fftw_destroy_plan(plan);
        fftw_cleanup();

	
        // TODO Calculating array of kx
       
    	for (int i = 0; i < (nc / 2 + 1); i++){
            kx[i] = i * delkx;
        }
        for (int i = -nc / 2; i < 0; i++){
            kx[i + nc + 1] = i * delkx;

            //For nyquist frequency
            if((nc + 1)%2 == 0){
                kx[nc/2 + 1] = nc * delkx / 2;
            }
        }
        
	// TODO Calculating Electric field in k space

  	for (int j = 0; j < nc + 1; j++){
            if(kx[j] != 0){
        	E_k[j][REAL] =  rho_k[j][IMAG] / kx[j];
        	E_k[j][IMAG] = -rho_k[j][REAL] / kx[j];
            }
            else{
            	E_k[j][REAL] = 0.0; //-rho_k[j][REAL]/kx[j+1];
            	E_k[j][IMAG] = 0.0; //-rho_k[j][IMAG]/kx[j+1];
            }
        }
	
	
	// TODO Calculating Electric field in real space
	
	fftw_plan plan_i = fftw_plan_dft_1d(nc+1, E_k, E_x, FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_i);

        //Do some cleaning
        fftw_destroy_plan(plan_i);
        fftw_cleanup();
       

	//Getting rid of unnormalisation due to fft
        for (int i=0; i<nc+1; i++){
        Ex_norm[i][REAL]=E_x[i][REAL]/(nc+1);
        Ex_norm[i][IMAG]=E_x[i][IMAG]/(nc+1);
        }
	
	//Storing real value of Electric field into E_field
	for(int j=0; j<nc+1; j++){
        E_field[j]=Ex_norm[j][REAL];
	}
 
     return 0;

}    

int half_velocity(){//Dev Harihar
    double Eph; //store interpolated electric field at the particle location
    //loop over particles
    for(int i = 0; i < ncols; i++){

        //TODO update velocity_x array to the half time step using the electric field
    	
	int j = (int)(ParArray[0][i] / dx);
    	Eph   = (E_field[j] * (xgrid[j + 1] - ParArray[0][i]) / dx) + (E_field[j + 1] * (ParArray[0][i] - xgrid[j]) / dx);
    	
	velocity_x[i] = v_0 + (dt / 2) * Eph * (ParArray[1][i] / ParArray[2][i]);

    }
    return 0;

}

int particle_push(){//dev HArihar
    //cout<<"Particle Move"<<endl;
    double Ep;

    //loop over particles
    for(int i = 0; i < ncols; i++){

        //TODO update particle position and velocity (which one should be done first?)

        int j = (int)(ParArray[0][i] / dx);     //getting grid index of particle

        //finding force on particle to get new velocity
        Ep = (E_field[j] * (xgrid[j + 1] - ParArray[0][i]) / dx) + (E_field[j + 1] * (ParArray[0][i] - xgrid[j]) / dx);

        ParArray[0][i] += dt * velocity_x[i];
        
	//new velocity
        velocity_x[i]  += dt * Ep * (ParArray[1][i] / ParArray[2][i]);
        

	//TODO apply periodic boundary conditions to particles which go out of the box

	if(ParArray[0][i] > Lx){
            ParArray[0][i] = ParArray[0][i] - Lx;
        }
        if(ParArray[0][i] < 0){
            ParArray[0][i] = ParArray[0][i] + Lx;
        }

    }
    return 0;
}
