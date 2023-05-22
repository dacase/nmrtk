/***************************************************************************************  
 *   myNamespaces.h, define the namespaces in this package.                            *
 *                                                                                     *  
 *                                                                                     *
 **************************************************************************************/
#ifndef MYNAMESPACE_H

#define MYNAMESPACE_H

namespace nmrspace
{
  const double Mass_C = 12.011;            // g/mol
  const double Mass_N = 14.0067;           // g/mol
  const double Mass_O = 15.9994;           // g/mol
  const double Mass_H = 1.008;             // g/mol
  const double Mass_S = 32.066;            // g/mol
  const double Mass_P = 30.97376;          // g/mol

  const double H_BAR  = 1.05457266e-34;    // Js
  const double MU_0   = 12.56e-7;          // Vs/(Am)
  const double PI     = 3.141592;          

  const double Gamma_C =  67.2274e6;       // rad/sT
  const double Gamma_N = -27.101968e6;     // rad/sT
  const double Gamma_O = -36.245648e6;     // rad/sT
  const double Gamma_H = 267.522128e6;     // rad/sT
  const double Gamma_P = 108.2358e6;       // rad/sT

  const double LFreq_C =  10.705e6;        // hz/T
  const double LFreq_N =   4.315e6;        // hz/T
  const double LFreq_O =   5.772e6;        // hz/T
  const double LFreq_H =  42.577e6;        // hz/T
  const double LFreq_P =  17.235e6;        // hz/T

  const double Reff_C  = 1.09e-10;         // m,  default value for C-H vector 
  const double Reff_N  = 1.01e-10;         // m,  default value for N-H vector 

  const double Delta_C =  25e-6;           
  const double Delta_N = 170e-6;             

}

#endif
