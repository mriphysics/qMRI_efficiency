#include <iostream>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include "mex.h"

const std::complex<double> i1(0.0, 1.0);

Eigen::Matrix3cd T;


void update_T(double a, double p)
{
    T(0,0) = pow(cos(a/2),2);
    T(0,1) = pow(sin(a/2),2)*(cos(2*p)+1.0*i1*sin(2*p));
    T(0,2) = sin(a)*(sin(p)-1.0*i1*cos(p));
    T(1,0) = std::conj(T(0,1));
    T(1,1) = T(0,0);
    T(1,2) = std::conj(T(0,2));
    T(2,0) = -0.5*sin(a)*(sin(p)+1.0*i1*cos(p));
    T(2,1) = std::conj(T(2,0));
    T(2,2) = cos(a);
}

void EPG_GRE(double *signal_real, double *signal_imag, double *path_real, double *path_imag, double *N, double *FA, double *RF, double *TR, double *T1, double *T2, double *M0)
{
    int n = (int) N[0];

    //max order of EPG simulation
    int KMAX = n;
    if(n>25){KMAX = 25;}

    ///// MEMORY ALLOCATION	
    // relaxation operators
    double E1 = exp(-TR[0]/T1[0]);
    double E2 = exp(-TR[0]/T2[0]);

    // signal vectors
    Eigen::VectorXcd F_plus(KMAX);
    F_plus.setZero();

    Eigen::VectorXcd F_minus(KMAX);
    F_minus.setZero();

    Eigen::VectorXcd Z(KMAX);
    Z.setZero();

    Eigen::Vector3cd aux_states(0, 0, 0);


    /////SIGNAL SIMULATION
    Z(0) = M0[0];
    
    // transition matrix
    update_T(FA[0], RF[0]);

    //apply first pulse
    aux_states(0) = F_plus(0);
    aux_states(1) = F_minus(0);
    aux_states(2) = Z(0);
    aux_states = T * aux_states;
    F_plus(0) = aux_states(0);
    F_minus(0) = aux_states(1);
    Z(0) = aux_states(2);

    //save F_minus_zero >> this is the measured signal
    signal_real[0] = F_minus(0).real();
    signal_imag[0] = F_minus(0).imag();

    //simulate rest of the sequence
    for (int j=1; j<n; j++){

	    //apply E
    	for (int k=0; k<(std::min(j,KMAX)); k++)
    	{
    		F_plus(k) = F_plus(k) * E2;
    		F_minus(k) = F_minus(k) * E2;
    		Z(k) = Z(k) * E1; 
    	}
    	Z(0) += (1.0 - E1)*M0[0];

    	//apply S
    	for (int k=(std::min(KMAX-1,j)); k>0; k--){ F_plus(k) = F_plus(k-1); }
    	F_plus(0) = std::conj(F_minus(1));
    	for (int k=0; k<(std::min(KMAX-1,j)); k++){ F_minus(k) = F_minus(k+1); }
    	F_minus(std::min(KMAX-1,j)) = 0.0 + 0.0*i1;
	
    	//apply T
        update_T(FA[j], RF[j]);
    	for (int k=0; k<(std::min(KMAX,j+1)); k++) 
    	{
    		aux_states(0) = F_plus(k);
    		aux_states(1) = F_minus(k);
    		aux_states(2) = Z(k);
    		aux_states = T * aux_states;
    
    		F_plus(k) = aux_states(0);
    		F_minus(k) = aux_states(1);
    		Z(k) = aux_states(2);
    	}

            
    	//save F_minus_zero
    	signal_real[j] = F_minus(0).real();
        signal_imag[j] = F_minus(0).imag();

        //save last states 
        if(j==(n-1))
        {
            for (int k=0; k<KMAX; k++)
            {
                path_real[3*k] = F_plus(k).real();
                path_real[3*k+1] = F_minus(k).real();
                path_real[3*k+2] = Z(k).real();
                path_imag[3*k] = F_plus(k).imag();
                path_imag[3*k+1] = F_minus(k).imag();
                path_imag[3*k+2] = Z(k).imag();
            }
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int rhs, const mxArray *prhs[])
{
    /* Declare arrays necessary to compute signals */ 
    double *N;
    double *FA; 
    double *RF; 
    double *TR; 
    double *T1;
    double *T2;
    double *M0;

    double *signal_real;
    double *signal_imag;
    double *path_real;
    double *path_imag;

    /* Check for proper number of arguments */
    
    /* Create pointers to the inputs */
    N = mxGetPr(prhs[0]);
    FA = mxGetPr(prhs[1]); 
    RF = mxGetPr(prhs[2]); 
    TR = mxGetPr(prhs[3]);
    T1 = mxGetPr(prhs[4]); 
    T2 = mxGetPr(prhs[5]);
    M0 = mxGetPr(prhs[6]);

    /* Create pointers to the outputs */
    plhs[0] = mxCreateDoubleMatrix((mwSize)*N, 1, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(3, 25, mxCOMPLEX);
    
    // Check this link to see how to return complex arrays:
    // http://matlab.izmiran.ru/help/techdoc/matlab_external/ch04cre9.html
    signal_real = mxGetPr(plhs[0]);
    signal_imag = mxGetPi(plhs[0]);
    path_real = mxGetPr(plhs[1]);
    path_imag = mxGetPi(plhs[1]);
        
    /* Call the computational routine */
    EPG_GRE(signal_real, signal_imag, path_real, path_imag, N, FA, RF, TR, T1, T2, M0);


}




