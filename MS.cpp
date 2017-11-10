/*----------------------------------------------------------------------------80
    Author: Wanshun DU <dws.en.france@gmail.com>
------------------------------------------------------------------------------*/

#include "math.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 8)
        mexErrMsgTxt("apply_pbc: Wrong number of input arguments.\n");
    if (nlhs > 8)   
        mexErrMsgTxt("apply_pbc: Too many output argumnents.\n");

    #define f_IN  prhs[0]
    #define fx_IN  prhs[1]
    #define fy_IN  prhs[2]
    #define fxx_IN  prhs[3]
    #define fyy_IN  prhs[4]
    #define fxy_IN  prhs[5]
    #define mu_IN  prhs[6]
    #define T_IN  prhs[7]
    #define chi_spin_OUT plhs[0]
    #define chi_LP_OUT plhs[1]
    #define chi_omega_a_OUT plhs[2]
    #define chi_omega_b_OUT plhs[3]
    #define chi_g_OUT plhs[4]
    #define chi_tildeg_a_OUT plhs[5]
    #define chi_tildeg_b_OUT plhs[6]
    #define chi_tildeg_c_OUT plhs[7]

    int N = mxGetN(mu_IN);

    chi_spin_OUT = mxCreateDoubleMatrix(1, N, mxREAL);
    chi_LP_OUT = mxCreateDoubleMatrix(1, N, mxREAL);
    chi_omega_a_OUT = mxCreateDoubleMatrix(1, N, mxREAL);
    chi_omega_b_OUT = mxCreateDoubleMatrix(1, N, mxREAL);
    chi_g_OUT = mxCreateDoubleMatrix(1, N, mxREAL);
    chi_tildeg_a_OUT = mxCreateDoubleMatrix(1, N, mxREAL);
    chi_tildeg_b_OUT = mxCreateDoubleMatrix(1, N, mxREAL);
    chi_tildeg_c_OUT = mxCreateDoubleMatrix(1, N, mxREAL);

    double *f = mxGetPr(f_IN);
    double *fx = mxGetPr(fx_IN);
    double *fy = mxGetPr(fy_IN);
    double *fxx = mxGetPr(fxx_IN);
    double *fyy = mxGetPr(fyy_IN);
    double *fxy = mxGetPr(fxy_IN);
    double *mu_I = mxGetPr(mu_IN);
    double *T_I = mxGetPr(T_IN);
        double T = T_I[0];
        double T_1 = 1 / T;

    double *chi_spin_O = mxGetPr(chi_spin_OUT);
    double *chi_LP_O = mxGetPr(chi_LP_OUT);
    double *chi_omega_a_O = mxGetPr(chi_omega_a_OUT);
    double *chi_omega_b_O = mxGetPr(chi_omega_b_OUT);
    double *chi_g_O = mxGetPr(chi_g_OUT);
    double *chi_tildeg_a_O = mxGetPr(chi_tildeg_a_OUT);
    double *chi_tildeg_b_O = mxGetPr(chi_tildeg_b_OUT);
    double *chi_tildeg_c_O = mxGetPr(chi_tildeg_c_OUT);

    double f_0 = 0;
    double fx_0 = 0;
    double fy_0 = 0;
    double fxx_0 = 0;
    double fyy_0 = 0;
    double fxy_0 = 0;

    double ee = sqrt ( f[0] * f[0] + f[1] * f[1] + 0.000000000001 );
    double ee_1 = 1 / ee;
    double ex = ( fx[0] * f[0] + fx[1] * f[1] ) * ee_1;
    double ey = ( fy[0] * f[0] + fy[1] * f[1] ) * ee_1;

    double nx[2], ny[2];
    for (int i = 0; i < 2; i++)
    {
        nx[i] = ( fx[i] - ex * f[i] * ee_1 ) * ee_1;
        ny[i] = ( fy[i] - ey * f[i] * ee_1 ) * ee_1;
    }

    double ep_x = fx_0 + ex;
    double em_x = fx_0 - ex;
    double ep_y = fy_0 + ey;
    double em_y = fy_0 - ey;

    double gxx = ( nx[0] * nx[0] + nx[1] * nx[1] ) / 4;
    double gyy = ( ny[0] * ny[0] + ny[1] * ny[1] ) / 4;
    double gxy = ( nx[0] * ny[0] + nx[1] * ny[1] ) / 4;

    double exx = ( fxx[0] * f[0] + fxx[1] * f[1] ) * ee_1 + 4 * ee * gxx;
    double eyy = ( fyy[0] * f[0] + fyy[1] * f[1] ) * ee_1 + 4 * ee * gyy;
    double exy = ( fxy[0] * f[0] + fxy[1] * f[1] ) * ee_1 + 4 * ee * gxy;

    double nxx[2], nyy[2], nxy[2];
    for (int i = 0; i < 2; i++)
    {
        nxx[i] = ( fxx[i] - 2 * ex * nx[i] - exx * f[i] * ee_1 ) * ee_1;
        nyy[i] = ( fyy[i] - 2 * ey * ny[i] - eyy * f[i] * ee_1 ) * ee_1;
        nxy[i] = ( fxy[i] - ex * ny[i] - ey * nx[i] - exy * f[i] * ee_1 ) * ee_1;
    }

    double gxxy = ( nxx[0] * ny[0] + nxx[1] * ny[1] ) / 4;
    double gyyx = ( nyy[0] * nx[0] + nyy[1] * nx[1] ) / 4;
    double gxyx = ( nxy[0] * nx[0] + nxy[1] * nx[1] ) / 4;
    double gxyy = ( nxy[0] * ny[0] + nxy[1] * ny[1] ) / 4;

    double gxxyy = ( nxx[0] * nyy[0] + nxx[1] * nyy[1] ) / 4;
    double gxyxy = ( nxy[0] * nxy[0] + nxy[1] * nxy[1] ) / 4;
            
    double exx_p = fxx_0 + exx;
    double eyy_p = fyy_0 + eyy;
    double exy_p = fxy_0 + exy;
    double exx_m = fxx_0 - exx;
    double eyy_m = fyy_0 - eyy;
    double exy_m = fxy_0 - exy;
           
    double LP_p = exx_p * eyy_p - exy_p * exy_p;
    double LP_m = exx_m * eyy_m - exy_m * exy_m;

    double V_omega = -ee * 4 * ( gxx * gyy - gxy * gxy );
    double V_g = -0.5 * ( ex * ( gyyx - gxyy ) + ey * ( gxxy - gxyx ) );
    double P_g = -0.5 * 2 * ee * ( gxxyy - gxyxy );
            
    double V_tildeg = ( fx_0 * fx_0 ) * gyy + ( fy_0 * fy_0 ) * gxx - 2 * fx_0 * fy_0 * gxy;
    double U_tildeg_a = ( fx_0 * ex ) * gyy + ( fy_0 * ey ) * gxx - fx_0 * ey *gxy - fy_0 * ex *gxy;
    double U_tildeg_b = fxx_0 * gyy + fyy_0 * gxx - 2 * fxy_0 * gxy;
    double P_tildeg = fx_0 * ( gyyx - gxyy ) + fy_0 * ( gxxy - gxyx );

    for (int j = 0; j < N; j++)
    {
    	double mu = mu_I[j];
    	    double ep = f_0 + ee - mu;
    	double em = f_0 - ee - mu;

    	double np_F = 1 / ( exp ( ep * T_1 ) + 1 ); 
    	double nm_F = 1 / ( exp ( em * T_1 ) + 1 ); 
            
    	double np_F1 = np_F * ( np_F - 1 ) * T_1;
    	double nm_F1 = nm_F * ( nm_F - 1 ) * T_1;
            
    chi_spin_O[j] = chi_spin_O[j] + np_F1 + nm_F1;
    chi_LP_O[j] = chi_LP_O[j] +( LP_m * nm_F1 + LP_p * np_F1 ) / 12;
    chi_omega_a_O[j] = chi_omega_a_O[j] + ( nm_F1 + np_F1 ) * ee * V_omega;
    chi_omega_b_O[j] = chi_omega_b_O[j] + ( nm_F - np_F ) * V_omega;
    chi_g_O[j] = chi_g_O[j] + ( nm_F - np_F ) * ( 2 * V_g + P_g );
    chi_tildeg_a_O[j] = chi_tildeg_a_O[j] + ( nm_F - np_F ) * ee_1 * V_tildeg;
    chi_tildeg_b_O[j] = chi_tildeg_b_O[j] - ( nm_F - np_F ) * ( P_tildeg + 2 * U_tildeg_b );
    chi_tildeg_c_O[j] = chi_tildeg_c_O[j] + ( nm_F1 + np_F1 ) * V_tildeg + ( np_F1 - nm_F1 ) * U_tildeg_a;
    }

}