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
void mex_isochromat_CRLB_theta(double *signal_real, double *signal_imag, double *dmdT1_real, double *dmdT1_imag, double *dmdT2_real, double *dmdT2_imag, double *Mz, double *rCRLB, double *N, double *theta, double *RF, double *TRssfp, double *T1, double *T2, double *M0, double *B0, double *phi0)
{

    double R1 = 1 / T1[0];
	double R2 = 1 / T2[0];

	int length = (int) N[0];

    ///// MEMORY ALLOCATION AND INITIALIZATION	
    // signal vectors
    Eigen::Vector4d Mplus(0, 0, 0, 1);  //Magnetization vector after the RF excitation
    Eigen::Vector4d Mminus(0, 0, 0, 1); //Magnetization vector before the RF excitation
    Mminus(2) = M0[0];

    Eigen::VectorXd aux_FA(length);
    aux_FA.setZero();
    Eigen::VectorXd TR(length);
    TR.setOnes();
    TR *= TRssfp[0];
    Eigen::VectorXd TE(length);
    TE.setOnes();
    TE *= TRssfp[0]/2.0;
    
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
    for (int j=0; j<length; j++){

        //transformation from theta to alpha
        if (j==0) {aux_FA(j) = theta[j];}
        else {aux_FA(j) = theta[j] + theta[j-1];}

        //apply rotations
        update_Rz(-RF[j]); update_Rx(aux_FA(j));
        Mplus = Rx * Rz * Mminus; 
        Mplus_dOdT1 = Rx * Rz * Mminus_dOdT1; 
        Mplus_dOdT2 = Rx * Rz * Mminus_dOdT2; 
        Mplus_dOdM0 = Rx * Rz * Mminus_dOdM0;
        Mplus_dOdB0 = Rx * Rz * Mminus_dOdB0;
        Mplus_dOdP0 = Rx * Rz * Mminus_dOdP0;
        update_Rz(RF[j]);
        Mplus = Rz * Mplus;
        Mplus_dOdT1 = Rz * Mplus_dOdT1; 
        Mplus_dOdT2 = Rz * Mplus_dOdT2; 
        Mplus_dOdM0 = Rz * Mplus_dOdM0;      
        Mplus_dOdB0 = Rz * Mplus_dOdB0;
        Mplus_dOdP0 = Rz * Mplus_dOdP0;

        //CRLB calculation
        FPangle = 2*PI*B0[0] * TE(j) * pow(10.0,-3.0);
        update_Rz(FPangle + phi0[0]);
        Mminus       = Rz * Mplus;
        Mminus_dOdT1 = Rz * Mplus_dOdT1; 
        Mminus_dOdT2 = Rz * Mplus_dOdT2; 
        Mminus_dOdM0 = Rz * Mplus_dOdM0;       
        
        update_dRz(FPangle + phi0[0], 2*PI*TE(j)*pow(10.0,-3.0) );       
        Mminus_dOdB0 = dRz * Mplus + Rz * Mplus_dOdB0;
        
        update_dRz(FPangle + phi0[0], 1.0);
        Mminus_dOdP0 = dRz * Mplus + Rz * Mplus_dOdP0;
        
        dmdT1 = exp(-R2 * TE(j)) * (Mminus_dOdT1(0) + i1*Mminus_dOdT1(1)); 
        dmdT2 = TE(j) * exp(-R2 * TE(j)) * pow(R2,2) * (Mminus(0) + i1*Mminus(1)) + exp(-R2 * TE(j)) * (Mminus_dOdT2(0) + i1*Mminus_dOdT2(1));
        dmdM0 = exp(-R2 * TE(j)) * (Mminus_dOdM0(0) + i1*Mminus_dOdM0(1));
        dmdB0 = exp(-R2 * TE(j)) * (Mminus_dOdB0(0) + i1*Mminus_dOdB0(1)); 
        dmdP0 = exp(-R2 * TE(j)) * (Mminus_dOdP0(0) + i1*Mminus_dOdP0(1));
        
        add2FIM(); 

        signal_real[j] = exp(-R2 * TE(j)) * Mminus(0);
        signal_imag[j] = exp(-R2 * TE(j)) * Mminus(1);
        dmdT1_real[j]  = dmdT1.real();
        dmdT1_imag[j]  = dmdT1.imag();
        dmdT2_real[j]  = dmdT2.real();
        dmdT2_imag[j]  = dmdT2.imag();
        Mz[j] = Mminus(2) * exp(-R1 * TR(j)) + (1 - exp(-R1 * TR(j))) * std::abs(M0[0]);

        //apply relaxation
        update_E(R1, R2, TR(j), std::abs(M0[0]));
    	update_dEdT1(R1, TR(j), std::abs(M0[0]));
        update_dEdT2(R2, TR(j));
        update_dEdM0(R1, TR(j));
        Mminus = E * Mplus;
        Mminus_dOdT1 = E * Mplus_dOdT1 + dEdT1 * Mplus;
        Mminus_dOdT2 = E * Mplus_dOdT2 + dEdT2 * Mplus;
        Mminus_dOdM0 = E * Mplus_dOdM0 + dEdM0 * Mplus;
        Mminus_dOdB0 = E * Mplus_dOdB0;
        Mminus_dOdP0 = E * Mplus_dOdP0;

        //apply rotation due to free precession
        FPangle = 2*PI*B0[0] * TR(j) * pow(10.0,-3.0);
        update_Rz(FPangle);
        update_dRz(FPangle, 2*PI*TR(j)*pow(10.0,-3.0) );
        Mminus = Rz * Mminus;
        Mminus_dOdT1 = Rz * Mminus_dOdT1; 
        Mminus_dOdT2 = Rz * Mminus_dOdT2; 
        Mminus_dOdM0 = Rz * Mminus_dOdM0;       
        Mminus_dOdB0 = dRz * E * Mplus + Rz * Mminus_dOdB0;
        Mminus_dOdP0 = Rz * Mminus_dOdP0;

}

    V = F.inverse();
    rCRLB[0] = V(0,0);
    rCRLB[1] = V(1,1);
    rCRLB[2] = V(2,2);
    rCRLB[3] = V(3,3);
    rCRLB[4] = V(4,4);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Declare arrays necessary to compute signals */ 
    double *N;         //units []
    double *theta;     //units [rad]
    double *RF;        //units [rad]
    double *TRssfp;    //units [ms]
    double *T1;        //units [ms]
    double *T2;        //units [ms]
    double *M0;        //a.u.
    double *B0;        //units [Hz]
    double *phi0;      //units [rad]

    double *signal_real;
    double *signal_imag;
    double *dmdT1_real;
    double *dmdT1_imag;   
    double *dmdT2_real;
    double *dmdT2_imag;     
    double *Mz;
    double *rCRLB;

    /* Check for proper number of arguments */
        
    /* Create pointers to the inputs */
    N      = mxGetPr(prhs[0]);
    theta  = mxGetPr(prhs[1]); 
    RF     = mxGetPr(prhs[2]); 
    TRssfp = mxGetPr(prhs[3]);
    T1     = mxGetPr(prhs[4]); 
    T2     = mxGetPr(prhs[5]);
    M0     = mxGetPr(prhs[6]);
    B0     = mxGetPr(prhs[7]);
    phi0   = mxGetPr(prhs[8]);

    /* Create pointers to the outputs */
    plhs[0] = mxCreateDoubleMatrix((mwSize)*N, 1, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix((mwSize)*N, 1, mxCOMPLEX);    
    plhs[2] = mxCreateDoubleMatrix((mwSize)*N, 1, mxCOMPLEX);
    plhs[3] = mxCreateDoubleMatrix((mwSize)*N, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, 5, mxREAL);
     
    // Check this link to see how to return complex arrays:
    // http://matlab.izmiran.ru/help/techdoc/matlab_external/ch04cre9.html
    signal_real = mxGetPr(plhs[0]);
    signal_imag = mxGetPi(plhs[0]);
    dmdT1_real  = mxGetPr(plhs[1]);
    dmdT1_imag  = mxGetPi(plhs[1]);
    dmdT2_real  = mxGetPr(plhs[2]);
    dmdT2_imag  = mxGetPi(plhs[2]);
    Mz          = mxGetPr(plhs[3]);
    rCRLB       = mxGetPr(plhs[4]);
       
    /* Call the computational routine */
    mex_isochromat_CRLB_theta(signal_real, signal_imag, dmdT1_real, dmdT1_imag, dmdT2_real, dmdT2_imag, Mz, rCRLB, N, theta, RF, TRssfp, T1, T2, M0, B0, phi0);

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





