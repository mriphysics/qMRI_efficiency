#include <iostream>
#include <complex>
#include <cmath>
#include <Eigen/Dense>
#include "mex.h"

const std::complex<double> i1(0.0, 1.0);
typedef Eigen::Matrix<double, 2, 3> Matrix2x3d;

Eigen::Matrix3cd T;
Eigen::Vector3cd aux_states(0, 0, 0);
Eigen::Matrix3d F;
Eigen::Matrix3d V;
Matrix2x3d J;
std::complex<double> dmdT1;
std::complex<double> dmdT2;
std::complex<double> dmdM0;

 

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

void add2FIM()
{
    J(0,0) = dmdT1.real();
    J(1,0) = dmdT1.imag();
    J(0,1) = dmdT2.real();
    J(1,1) = dmdT2.imag(); 
    J(0,2) = dmdM0.real();
    J(1,2) = dmdM0.imag();  
    F     += J.transpose() * J;
}

void EPG_GRE_efficiency(double *eff, double *signal_real, double *signal_imag, double *N, double *FA, double *RF, double *TR, double *TE, double *T1, double *T2)
{
    int n = (int) N[0];

    //max order of EPG simulation
    int KMAX = n;
    if(n>25){KMAX = 25;}

    ///// MEMORY ALLOCATION	
    double E1;
    double E2;
    double E2_TE;
    double dEdT1;
    double dEdT2;
    double dEdT2_TE;
    double Tacq = 0.0;

    // fisher information matrix
    F.setZero();

    // signal vectors
    Eigen::VectorXcd Fplus(KMAX);
    Fplus.setZero();
    Eigen::VectorXcd Fminus(KMAX);
    Fminus.setZero();
    Eigen::VectorXcd Z(KMAX);
    Z.setZero();
    Eigen::VectorXcd prev_Fplus(KMAX);
    prev_Fplus.setZero();
    Eigen::VectorXcd prev_Fminus(KMAX);
    prev_Fminus.setZero();
    Eigen::VectorXcd prev_Z(KMAX);
    prev_Z.setZero();

    // for CRLB calculation
    Eigen::VectorXcd Fplus_dOdT1(KMAX);
    Fplus_dOdT1.setZero();
    Eigen::VectorXcd Fminus_dOdT1(KMAX);
    Fminus_dOdT1.setZero();   
    Eigen::VectorXcd Z_dOdT1(KMAX);
    Z_dOdT1.setZero();
    Eigen::VectorXcd prev_Fplus_dOdT1(KMAX);
    prev_Fplus_dOdT1.setZero();
    Eigen::VectorXcd prev_Fminus_dOdT1(KMAX);
    prev_Fminus_dOdT1.setZero();
    Eigen::VectorXcd prev_Z_dOdT1(KMAX);
    prev_Z_dOdT1.setZero();
    
    Eigen::VectorXcd Fplus_dOdT2(KMAX);
    Fplus_dOdT2.setZero();
    Eigen::VectorXcd Fminus_dOdT2(KMAX);
    Fminus_dOdT2.setZero();   
    Eigen::VectorXcd Z_dOdT2(KMAX);
    Z_dOdT2.setZero();
    Eigen::VectorXcd prev_Fplus_dOdT2(KMAX);
    prev_Fplus_dOdT2.setZero();
    Eigen::VectorXcd prev_Fminus_dOdT2(KMAX);
    prev_Fminus_dOdT2.setZero();
    Eigen::VectorXcd prev_Z_dOdT2(KMAX);
    prev_Z_dOdT2.setZero();
    
    Eigen::VectorXcd Fplus_dOdM0(KMAX);
    Fplus_dOdM0.setZero();
    Eigen::VectorXcd Fminus_dOdM0(KMAX);
    Fminus_dOdM0.setZero();   
    Eigen::VectorXcd Z_dOdM0(KMAX);
    Z_dOdM0.setZero();
    Eigen::VectorXcd prev_Fplus_dOdM0(KMAX);
    prev_Fplus_dOdM0.setZero();
    Eigen::VectorXcd prev_Fminus_dOdM0(KMAX);
    prev_Fminus_dOdM0.setZero();
    Eigen::VectorXcd prev_Z_dOdM0(KMAX);
    prev_Z_dOdM0.setZero();       

    /////CRLB SIMULATION
    Z(0)       = 1.0 + 0.0*i1;
    Z_dOdM0(0) = 1.0 + 0.0*i1;
    
    //simulate the sequence
    for (int j=0; j<n; j++){

        E1       = exp(-TR[j] / T1[0]);
        E2       = exp(-TR[j] / T2[0]);
        E2_TE    = exp(-TE[j] / T2[0]);
        dEdT1    = TR[j] * E1 / pow(T1[0],2);
        dEdT2    = TR[j] * E2 / pow(T2[0],2);
        dEdT2_TE = TE[j] * E2_TE / pow(T2[0],2);
        Tacq    += TR[j]*1e-3;

        //apply T
        update_T(FA[j], RF[j]);
    	for (int k=0; k<(std::min(KMAX,j+1)); k++) 
    	{
    		aux_states(0) = Fplus(k);
    		aux_states(1) = Fminus(k);
    		aux_states(2) = Z(k);
    		aux_states    = T * aux_states;
            
    		prev_Fplus(k)  = aux_states(0);
    		prev_Fminus(k) = aux_states(1);
    		prev_Z(k)      = aux_states(2);
            
            //T1
            aux_states(0) = Fplus_dOdT1(k);
  		    aux_states(1) = Fminus_dOdT1(k);
   		    aux_states(2) = Z_dOdT1(k);
   		    aux_states    = T * aux_states;
            
   		    prev_Fplus_dOdT1(k)  = aux_states(0);
   		    prev_Fminus_dOdT1(k) = aux_states(1);
   		    prev_Z_dOdT1(k)      = aux_states(2);    

            //T2
            aux_states(0) = Fplus_dOdT2(k);
  		    aux_states(1) = Fminus_dOdT2(k);
   		    aux_states(2) = Z_dOdT2(k);
   		    aux_states    = T * aux_states;
            
   		    prev_Fplus_dOdT2(k)  = aux_states(0);
   		    prev_Fminus_dOdT2(k) = aux_states(1);
   		    prev_Z_dOdT2(k)      = aux_states(2);

            //M0
            aux_states(0) = Fplus_dOdM0(k);
  		    aux_states(1) = Fminus_dOdM0(k);
   		    aux_states(2) = Z_dOdM0(k);
   		    aux_states    = T * aux_states;
            
   		    prev_Fplus_dOdM0(k)  = aux_states(0);
   		    prev_Fminus_dOdM0(k) = aux_states(1);
   		    prev_Z_dOdM0(k)      = aux_states(2);    
    	}

        signal_real[j] = prev_Fminus(0).real();
        signal_imag[j] = prev_Fminus(0).imag();

        //CRLB calculation
        dmdT1 = E2_TE * prev_Fminus_dOdT1(0); 
        dmdT2 = dEdT2_TE * prev_Fminus(0) + E2_TE * prev_Fminus_dOdT2(0);
        dmdM0 = E2_TE * prev_Fminus_dOdM0(0);
        add2FIM();

	    //apply E
    	for (int k=0; k<(std::min(j+1,KMAX)); k++)
    	{
    		Fplus(k)  = prev_Fplus(k) * E2;
    		Fminus(k) = prev_Fminus(k) * E2;
    		Z(k)      = prev_Z(k) * E1; 
            
            Fplus_dOdT1(k)  = prev_Fplus_dOdT1(k) * E2;
            Fminus_dOdT1(k) = prev_Fminus_dOdT1(k) * E2;
            Z_dOdT1(k)      = prev_Z(k) * dEdT1 + prev_Z_dOdT1(k) * E1;

            Fplus_dOdT2(k)  = prev_Fplus_dOdT2(k) * E2 + prev_Fplus(k) * dEdT2;
            Fminus_dOdT2(k) = prev_Fminus_dOdT2(k) * E2 + prev_Fminus(k) * dEdT2;
            Z_dOdT2(k)      = prev_Z_dOdT2(k) * E1;

            Fplus_dOdM0(k)  = prev_Fplus_dOdM0(k) * E2;
            Fminus_dOdM0(k) = prev_Fminus_dOdM0(k) * E2;
            Z_dOdM0(k)      = prev_Z_dOdM0(k) * E1;
       	}
    	Z(0)       += 1.0 - E1;
        Z_dOdT1(0) -= dEdT1;
        Z_dOdM0(0) += 1.0 - E1;

    	//apply S
    	for (int k=(std::min(KMAX-1,j+1)); k>0; k--)
        {
            Fplus(k)       = Fplus(k-1); 
            Fplus_dOdT1(k) = Fplus_dOdT1(k-1); 
            Fplus_dOdT2(k) = Fplus_dOdT2(k-1);        
            Fplus_dOdM0(k) = Fplus_dOdM0(k-1);
        }
    	Fplus(0)       = std::conj(Fminus(1));
        Fplus_dOdT1(0) = std::conj(Fminus_dOdT1(1));
        Fplus_dOdT2(0) = std::conj(Fminus_dOdT2(1));        
        Fplus_dOdM0(0) = std::conj(Fminus_dOdM0(1));
    	for (int k=0; k<(std::min(KMAX-1,j)); k++)
        { 
            Fminus(k)       = Fminus(k+1); 
            Fminus_dOdT1(k) = Fminus_dOdT1(k+1);           
            Fminus_dOdT2(k) = Fminus_dOdT2(k+1);
            Fminus_dOdM0(k) = Fminus_dOdM0(k+1);

        }
    	Fminus(std::min(KMAX-1,j))       = 0.0 + 0.0*i1;
		Fminus_dOdT1(std::min(KMAX-1,j)) = 0.0 + 0.0*i1;
        Fminus_dOdT2(std::min(KMAX-1,j)) = 0.0 + 0.0*i1;
		Fminus_dOdM0(std::min(KMAX-1,j)) = 0.0 + 0.0*i1;            
    }

    V = F.inverse();

    eff[0] = T1[0] / pow(V(0,0) * Tacq, 0.5); 
    eff[1] = T2[0] / pow(V(1,1) * Tacq, 0.5); 
 
}

void mexFunction(int nlhs, mxArray *plhs[], int rhs, const mxArray *prhs[])
{
    /* Declare arrays necessary to compute signals */ 
    double *N;
    double *FA; 
    double *RF; 
    double *TR; 
    double *TE;
    double *T1;
    double *T2;

    double *eff;
    double *signal_real;
    double *signal_imag;

    /* Check for proper number of arguments */
    
    /* Create pointers to the inputs */
    N  = mxGetPr(prhs[0]);
    FA = mxGetPr(prhs[1]); 
    RF = mxGetPr(prhs[2]); 
    TR = mxGetPr(prhs[3]);
    TE = mxGetPr(prhs[4]);
    T1 = mxGetPr(prhs[5]); 
    T2 = mxGetPr(prhs[6]);

    /* Create pointers to the outputs */
    plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)*N, 1, mxCOMPLEX);
    
    // Check this link to see how to return complex arrays:
    // http://matlab.izmiran.ru/help/techdoc/matlab_external/ch04cre9.html
    eff         = mxGetPr(plhs[0]);
    signal_real = mxGetPr(plhs[1]);
    signal_imag = mxGetPi(plhs[1]);
        
    /* Call the computational routine */
    EPG_GRE_efficiency(eff, signal_real, signal_imag, N, FA, RF, TR, TE, T1, T2);


}




