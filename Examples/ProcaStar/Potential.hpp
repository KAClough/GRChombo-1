/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "ADMFixedBGVars.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double mass;
        double self_interaction;
    };

    // class params
    params_t m_params;

    // add alias for metric vars
    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}


    //! Set the potential function for the proca field here
    template <class data_t, template <typename> class vars_t>
    void compute_stress_energy(data_t &rho_potential, 
			   Tensor<1, data_t> &Si_potential,
			   Tensor<2, data_t> &Sij_potential,
                           const vars_t<data_t> &vars,
                           const vars_t<Tensor<1, data_t>> &d1,
			   const Tensor<2, data_t> &gamma_UU) const
    {
	// defining some useful variables to ease the load 	

    	const double msquared = pow(m_params.mass, 2.0);
        const double c4 = m_params.self_interaction;

	data_t Asquared = 0;	
        const data_t Avec0_Resquared = vars.Avec0_Re * vars.Avec0_Re ;     	
    	FOR2(i, j)
    	{
        	Asquared += gamma_UU[i][j] * vars.chi * vars.Avec_Re[i] * vars.Avec[j];
    	}

	
	// massive terms 
	rho_potential = 0.5 * msquared * Avec0_Resquared 
		      + 0.5 * msquared * Asquared; 
	// Self interacting terms 
	rho_potential +=     c4 * msquared * Asquared * Asquared 
		       + 2. * c4 * msquared * Asquared * Avec0_Resquared 
		       - 3. * c4 * msquared * Avec0_Resquared * Avec0_Resquared;

    	FOR1(i)
    	{
				// massive terms 
        	Si_potential[i] = msquared * vars.Avec0_Re * vars.Avec_Re[i]
				// Self interacting terms 
				+ 4. * c4 * msquared * vars.Avec0_Re * vars.Avec_Re[i] * Asquared 
				- 4. * c4 * msquared * vars.Avec0_Re * vars.Avec_Re[i] * Avec0_Resquared;
    	}

	
    	FOR2(i, j)
    	{
		// massive terms 
        	Sij_potential[i][j] =
            	msquared * (vars.Avec_Re[i] * vars.Avec[j] +
                        0.5 * vars.h[i][j] / vars.chi * vars.Avec0_Re * vars.Avec0_Re);

            	Sij_potential[i][j] += - 0.5 * vars.h[i][j] / vars.chi *
                                 msquared * Asquared;
        	
    		// Self interacting terms 
        	Sij_potential[i][j] +=
            		  4. * c4 * msquared * vars.Avec_Re[i] * vars.Avec[j] * Asquared
		 	-      c4 * msquared * vars.h[i][j] / vars.chi * Asquared * Asquared
			- 4. * c4 * msquared * vars.Avec_Re[i] * vars.Avec[j] * Avec0_Resquared
			+ 2. * c4 * msquared * vars.h[i][j] / vars.chi * Asquared * Avec0_Resquared 
			-      c4 * msquared * vars.h[i][j] / vars.chi * Avec0_Resquared * Avec0_Resquared; 
	 	}

	
    }

    //! Set the potential function for the proca field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &dVdA_Re, data_t &dVdA_Im,
			   data_t &dAvec0dt_Re, data_t &dAvec0dt_Im,
                           const vars_t<data_t> &vars,
                           const vars_t<Tensor<1, data_t>> &d1) const
    {
        // calculate full spatial christoffel symbols and gamma^ij
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse(vars.h);
        const auto chris = compute_christoffel(d1.h, gamma_UU);
        const auto chris_phys = compute_phys_chris(d1.chi, vars.chi,
						   vars.h, gamma_UU,
						   chris.ULL);

        // for ease of reading
        const double c4 = m_params.self_interaction;
	
	Tensor<2, data_t> K_tensor;

	FOR2(i, j)
	{
	    K_tensor[i][j] = 1. / vars.chi * (vars.A[i][j] + 1. / 3. *
						vars.K * vars.h[i][j]);
	}

        // Here we are defining often used terms
        // DA[i][j] = D_i A_j
        Tensor<2, data_t> DA;
        FOR2(i, j)
        {
            DA[i][j] = d1.Avec[j][i];
            FOR1(k) { DA[i][j] += -chris_phys[k][i][j] * vars.Avec[k]; }
        }

        // DAScalar = D_i A^i
        data_t DA_scalar;
        DA_scalar = 0;
        FOR2(i, j) { DA_scalar += DA[i][j] * gamma_UU[i][j] * vars.chi; }

        // Xsquared = X^/mu X_/mu
        data_t Xsquared;
        Xsquared = -vars.Avec0_Re * vars.Avec0_Re;
        FOR2(i, j) { Xsquared += gamma_UU[i][j] * vars.chi * vars.Avec[j] * vars.Avec_Re[i]; }

        // C = 1 + 4 c4 A^k A_k - 12 c4 Avec0_Re^2
        data_t C = 1.0 - 12.0 * c4 * vars.Avec0_Re * vars.Avec0_Re;
        FOR2(i, j)
        {
            C += 4.0 * c4 * gamma_UU[i][j] * vars.chi * vars.Avec[j] * vars.Avec_Re[i];
        }

        // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
        dVdA_Re = pow(m_params.mass, 2.0) * (1.0 + 4.0 * c4 * Xsquared);

        // dphidt - for now the whole thing is here since it depends mainly
        // on the form of the potential - except the advection term which is in
        // the ProcaField code
        dAvec0dt_Re = 0;
        FOR2(i, j)
        {
            dAvec0dt_Re += -gamma_UU[i][j] * vars.chi * vars.Avec_Re[i] * d1.lapse[j];
        }
        // QUESTION: Should this be lapse * Z / C  or lapse * Z??
        dAvec0dt_Re += -vars.lapse * vars.Z / C;
        FOR4(i, j, k, l)
        {
            dAvec0dt_Re += -8.0 * c4 * vars.lapse / C * gamma_UU[i][k] *
                      gamma_UU[j][l] * vars.Avec_Re[i] * vars.Avec[j] * DA[k][l];
        }
        dAvec0dt_Re += vars.lapse / C * (1.0 + 4.0 * c4 * Xsquared) *
                  (vars.K * vars.Avec0_Re - DA_scalar);
        FOR1(i)
        {
            dAvec0dt_Re += 8.0 * c4 * vars.Avec0_Re * vars.lapse / C *
                      (vars.Evec[i] * vars.Avec_Re[i]);
        }
        FOR4(i, j, k, l)
        {
            dAvec0dt_Re += 8.0 * c4 * vars.Avec0_Re * vars.lapse / C *
                      (-K_tensor[i][j] * vars.Avec[k] *
                       vars.Avec[l] * gamma_UU[i][k] * gamma_UU[j][l]);
        }
        FOR2(i, j)
        {
            dAvec0dt_Re += 8.0 * c4 * vars.Avec0_Re * vars.lapse / C *
                      (2.0 * vars.Avec_Re[i] * d1.Avec0_Re[j] * gamma_UU[i][j] * vars.chi);
        }
    }
};

#endif /* POTENTIAL_HPP_ */
