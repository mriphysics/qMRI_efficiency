#include <iostream>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include "mex.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

typedef Eigen::Matrix<double, 2, 5> Matrix2x5d;
typedef Eigen::Matrix<double, 5, 5> Matrix5d;

Matrix5d F;
Matrix5d V;
Matrix2x5d J;

std::complex<double> dmdT1;
std::complex<double> dmdT2;
std::complex<double> dmdM0;
std::complex<double> dmdB0;
std::complex<double> dmdP0;
const std::complex<double> i1(0.0, 1.0);
const double PI = 3.14159265358979323846;
Eigen::Matrix4d Rz; 
Eigen::Matrix4d dRz;
Eigen::Matrix4d Rx; 
Eigen::Matrix4d E; 
Eigen::Matrix4d dEdT1; 
Eigen::Matrix4d dEdT2; 
Eigen::Matrix4d dEdM0;


void update_Rz(double p);

void update_dRz(double p, double dp);

void update_Rx(double a);

void update_E(double R1, double R2, double TR, double M0);

void update_dEdT1(double R1, double TR, double M0);

void update_dEdT2(double R2, double TR);

void update_dEdM0(double R1, double TR);

void add2FIM();


//Main Function that generates the EPG states given the acquisition parameters
void DE_isochromat_efficiency(double *f, double *N, double *theta, double *RF, double *TR0, double *T1, double *T2, double *M0, double *B0, double *P0)
{

    int length = (int) N[0];
    double R1 {1/T1[0]};
	double R2 {1/T2[0]};

    double Tacq = TR0[0] * N[0];
    int Nloops = (int) (10.0 * T1[0] / Tacq);
    if (Nloops<3) {Nloops = 3;}


    ///// MEMORY ALLOCATION AND INITIALIZATION	
    // signal vectors
    Eigen::Vector4d Mplus(0, 0, 0, 1);  //Magnetization vector after the RF excitation
    Eigen::Vector4d Mminus(0, 0, 0, 1); //Magnetization vector before the RF excitation
    Mminus(2) = M0[0];

    Eigen::VectorXd aux_FA(length);
    aux_FA.setZero();
    Eigen::VectorXd TR(length);
    TR.setOnes();
    TR *= TR0[0];
    Eigen::VectorXd TE(length);
    TE.setOnes();
    TE *= TR0[0]/2.0;
    
    // rotation and relaxation matrices
    Rz.setZero(); 
    Rz(2,2) = 1;
    Rz(3,3) = 1;
    
    Rx.setZero(); 
    Rx(0,0) = 1;
    Rx(3,3) = 1;
    
    E.setZero();
    E(3,3) = 1;
    dEdT1.setZero();
    dEdT2.setZero();
    dEdM0.setZero();
    
    // for CRLB calculation
    F.setZero();
    
    Eigen::Vector4d Mplus_dOdT1(0, 0, 0, 0);
    Eigen::Vector4d Mminus_dOdT1(0, 0, 0, 0);
    Eigen::Vector4d Mplus_dOdT2(0, 0, 0, 0);
    Eigen::Vector4d Mminus_dOdT2(0, 0, 0, 0);
    Eigen::Vector4d Mplus_dOdM0(0, 0, 0, 0);
    Eigen::Vector4d Mminus_dOdM0(0, 0, 0, 0);
    Mminus_dOdM0(2) = M0[0];
    Eigen::Vector4d Mplus_dOdB0(0, 0, 0, 0);
    Eigen::Vector4d Mminus_dOdB0(0, 0, 0, 0);
    Eigen::Vector4d Mplus_dOdP0(0, 0, 0, 0);
    Eigen::Vector4d Mminus_dOdP0(0, 0, 0, 0);  
 
    double FPangle;
         
    /////ISOCHROMAT DYNAMICS SIMULATION
    for (int j=0; j <(length*Nloops); j++){

        //transformation from theta to alpha
        if (j==0) {aux_FA(j) = theta[j];}
        else {aux_FA(j%length) = theta[j%length] + theta[(j-1)%length];}
 
        //Apply inversion: Mz -> -Mz , and spoiling: Mxy -> 0
        if(j%length==0){
            Mminus(0) = 0.0;
            Mminus(1) = 0.0;           
            Mminus(2) = -Mminus(2);

            Mminus_dOdT1(0) = 0.0;
            Mminus_dOdT1(1) = 0.0;
            Mminus_dOdT1(2) = -Mminus_dOdT1(2);
            Mminus_dOdT2(0) = 0.0;
            Mminus_dOdT2(1) = 0.0;
            Mminus_dOdT2(2) = -Mminus_dOdT2(2);           
            Mminus_dOdM0(0) = 0.0;
            Mminus_dOdM0(1) = 0.0;
            Mminus_dOdM0(2) = -Mminus_dOdM0(2);           
            Mminus_dOdB0(0) = 0.0;
            Mminus_dOdB0(1) = 0.0;
            Mminus_dOdB0(2) = -Mminus_dOdB0(2);           
            Mminus_dOdP0(0) = 0.0;
            Mminus_dOdP0(1) = 0.0;
            Mminus_dOdP0(2) = -Mminus_dOdP0(2);           
        }


        //apply rotations
        update_Rz(-RF[j%length]); update_Rx(aux_FA(j%length));
        Mplus = Rx * Rz * Mminus; 
        Mplus_dOdT1 = Rx * Rz * Mminus_dOdT1; 
        Mplus_dOdT2 = Rx * Rz * Mminus_dOdT2; 
        Mplus_dOdM0 = Rx * Rz * Mminus_dOdM0;
        Mplus_dOdB0 = Rx * Rz * Mminus_dOdB0;
        Mplus_dOdP0 = Rx * Rz * Mminus_dOdP0;
        update_Rz(RF[j%length]);
        Mplus = Rz * Mplus;
        Mplus_dOdT1 = Rz * Mplus_dOdT1; 
        Mplus_dOdT2 = Rz * Mplus_dOdT2; 
        Mplus_dOdM0 = Rz * Mplus_dOdM0;      
        Mplus_dOdB0 = Rz * Mplus_dOdB0;
        Mplus_dOdP0 = Rz * Mplus_dOdP0;

        //CRLB calculation
        if (j>=(Nloops-1)*length){
            FPangle = 2*PI*B0[0] * TE(j%length) * pow(10.0,-3.0);
            update_Rz(FPangle + P0[0]);
            Mminus       = Rz * Mplus;
            Mminus_dOdT1 = Rz * Mplus_dOdT1; 
            Mminus_dOdT2 = Rz * Mplus_dOdT2; 
            Mminus_dOdM0 = Rz * Mplus_dOdM0;       
             
            update_dRz(FPangle + P0[0], 2*PI*TE(j%length)*pow(10.0,-3.0) );       
            Mminus_dOdB0 = dRz * Mplus + Rz * Mplus_dOdB0;
        
            update_dRz(FPangle + P0[0], 1.0);
            Mminus_dOdP0 = dRz * Mplus + Rz * Mplus_dOdP0;
        
            dmdT1 = exp(-R2 * TE(j%length)) * (Mminus_dOdT1(0) + i1*Mminus_dOdT1(1)); 
            dmdT2 = TE(j%length) * exp(-R2 * TE(j%length)) * pow(R2,2) * (Mminus(0) + i1*Mminus(1)) + exp(-R2 * TE(j%length)) * (Mminus_dOdT2(0) + i1*Mminus_dOdT2(1));
            dmdM0 = exp(-R2 * TE(j%length)) * (Mminus_dOdM0(0) + i1*Mminus_dOdM0(1));
            dmdB0 = exp(-R2 * TE(j%length)) * (Mminus_dOdB0(0) + i1*Mminus_dOdB0(1)); 
            dmdP0 = exp(-R2 * TE(j%length)) * (Mminus_dOdP0(0) + i1*Mminus_dOdP0(1));
        
            add2FIM();
        }

        //apply relaxation
        update_E(R1, R2, TR(j%length), std::abs(M0[0]));
    	update_dEdT1(R1, TR(j%length), std::abs(M0[0]));
        update_dEdT2(R2, TR(j%length));
        update_dEdM0(R1, TR(j%length));
        Mminus = E * Mplus;
        Mminus_dOdT1 = E * Mplus_dOdT1 + dEdT1 * Mplus;
        Mminus_dOdT2 = E * Mplus_dOdT2 + dEdT2 * Mplus;
        Mminus_dOdM0 = E * Mplus_dOdM0 + dEdM0 * Mplus;
        Mminus_dOdB0 = E * Mplus_dOdB0;
        Mminus_dOdP0 = E * Mplus_dOdP0;

        //apply rotation due to free precession
        FPangle = 2*PI*B0[0] * TR(j%length) * pow(10.0,-3.0);
        update_Rz(FPangle);
        update_dRz(FPangle, 2*PI*TR(j%length)*pow(10.0,-3.0) );
        Mminus = Rz * Mminus;
        Mminus_dOdT1 = Rz * Mminus_dOdT1; 
        Mminus_dOdT2 = Rz * Mminus_dOdT2; 
        Mminus_dOdM0 = Rz * Mminus_dOdM0;       
        Mminus_dOdB0 = dRz * E * Mplus + Rz * Mminus_dOdB0;
        Mminus_dOdP0 = Rz * Mminus_dOdP0;
        
    }

    V = F.inverse();

    f[0] =  (V(0,0)*R1*R1  +  V(1,1)*R2*R2) * (Tacq*1e-3);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Declare arrays necessary to compute signals */ 
    double *N;         //units []
    double *theta;     //units [rad]
    double *RF;        //units [rad]
    double *TR;        //units [ms]
    double *T1;        //units [ms]
    double *T2;        //units [ms]
    double *M0;        //a.u.
    double *B0;        //units [Hz]
    double *P0;        //units [rad]

    double *f;

    /* Check for proper number of arguments */
        
    /* Create pointers to the inputs */
    N      = mxGetPr(prhs[0]);
    theta  = mxGetPr(prhs[1]); 
    RF     = mxGetPr(prhs[2]); 
    TR     = mxGetPr(prhs[3]);
    T1     = mxGetPr(prhs[4]); 
    T2     = mxGetPr(prhs[5]);
    M0     = mxGetPr(prhs[6]);
    B0     = mxGetPr(prhs[7]);
    P0     = mxGetPr(prhs[8]);

    /* Create pointers to the outputs */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
     
    // Check this link to see how to return complex arrays:
    // http://matlab.izmiran.ru/help/techdoc/matlab_external/ch04cre9.html
    f = mxGetPr(plhs[0]);
       
    /* Call the computational routine */
    DE_isochromat_efficiency(f, N, theta, RF, TR, T1, T2, M0, B0, P0);

}

/////AUXILIAR FUNCTIONS
//Update augmented rotation matrix around z-axis
void update_Rz(double p)
{
    Rz(0,0) = cos(p);
    Rz(0,1) = -sin(p);
    Rz(1,0) = sin(p);
    Rz(1,1) = cos(p);
}

//Update derivative of the augmented rotation matrix around z-axis
void update_dRz(double p, double dp)
{
    dRz(0,0) = -dp * sin(p);
    dRz(0,1) = -dp * cos(p);
    dRz(1,0) =  dp * cos(p);
    dRz(1,1) = -dp * sin(p);
}

//Update augmented rotation matrix around x-axis
void update_Rx(double a)
{
    Rx(1,1) = cos(a);
    Rx(1,2) = sin(a);
    Rx(2,1) = -sin(a);
    Rx(2,2) = cos(a);
}

//Update augmented relaxation matrix and its derivatives
void update_E(double R1, double R2, double TR, double M0)
{
	E(0,0) = exp(-R2 * TR);
	E(1,1) = exp(-R2 * TR);
	E(2,2) = exp(-R1 * TR);
	E(2,3) = (1 - exp(-R1 * TR)) * M0;
}

void update_dEdT1(double R1, double TR, double M0)
{
    dEdT1(2,2) = TR * exp(-R1 * TR) * pow(R1,2.0);
    dEdT1(2,3) = -dEdT1(2,2) * M0;
}

void update_dEdT2(double R2, double TR)
{
    dEdT2(0,0) = TR * exp(-R2 * TR) * pow(R2,2.0);
    dEdT2(1,1) = dEdT2(0,0);
}

void update_dEdM0(double R1, double TR)
{
    dEdM0(2,3) = 1 - exp(-R1 * TR);
}

//Update Fisher Information Matrix
void add2FIM()
{
    J(0,0) = dmdT1.real();
    J(1,0) = dmdT1.imag();
    J(0,1) = dmdT2.real();
    J(1,1) = dmdT2.imag(); 
    J(0,2) = dmdM0.real();
    J(1,2) = dmdM0.imag();  
    J(0,3) = dmdB0.real();
    J(1,3) = dmdB0.imag();
    J(0,4) = dmdP0.real();
    J(1,4) = dmdP0.imag();
    F += J.transpose() * J;
}





