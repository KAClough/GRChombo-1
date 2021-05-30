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
                           const vars_t<Tensor<1, data_t>> &d1) const
    {
	// defining some useful variables to ease the load 	

        using namespace TensorAlgebra;
    	const double msquared = pow(m_params.mass, 2.0);
        const double c4 = m_params.self_interaction;
        const auto gamma_UU = compute_inverse(vars.h);

	rho_potential = 0;	

	FOR1(i){Si_potential[i] = 0;}
	FOR2(i, j){Sij_potential[i][j] = 0;}

	data_t Asquared_Re = 0;	
        const data_t Avec0_Resquared = vars.Avec0_Re * vars.Avec0_Re ;     	
    	FOR2(i, j)
    	{
        	Asquared_Re += gamma_UU[i][j] * vars.chi * vars.Avec_Re[i] * vars.Avec_Re[j];
    	}

	
	// massive terms 
	rho_potential += 0.5 * msquared * Avec0_Resquared 
		      + 0.5 * msquared * Asquared_Re; 
	// Self interacting terms 
	rho_potential +=     c4 * msquared * Asquared_Re * Asquared_Re 
		       + 2. * c4 * msquared * Asquared_Re * Avec0_Resquared 
		       - 3. * c4 * msquared * Avec0_Resquared * Avec0_Resquared;

    	FOR1(i)
    	{
				// massive terms 
        	Si_potential[i] += msquared * vars.Avec0_Re * vars.Avec_Re[i]
				// Self interacting terms 
				+ 4. * c4 * msquared * vars.Avec0_Re * vars.Avec_Re[i] * Asquared_Re 
				- 4. * c4 * msquared * vars.Avec0_Re * vars.Avec_Re[i] * Avec0_Resquared;
    	}

	
    	FOR2(i, j)
    	{
		// massive terms 
        	Sij_potential[i][j] +=
            	msquared * (vars.Avec_Re[i] * vars.Avec_Re[j] +
                        0.5 * vars.h[i][j] / vars.chi * vars.Avec0_Re * vars.Avec0_Re);

            	Sij_potential[i][j] += - 0.5 * vars.h[i][j] / vars.chi *
                                 msquared * Asquared_Re;
        	
    		// Self interacting terms 
        	Sij_potential[i][j] +=
            		  4. * c4 * msquared * vars.Avec_Re[i] * vars.Avec_Re[j] * Asquared_Re
		 	-      c4 * msquared * vars.h[i][j] / vars.chi * Asquared_Re * Asquared_Re
			- 4. * c4 * msquared * vars.Avec_Re[i] * vars.Avec_Re[j] * Avec0_Resquared
			+ 2. * c4 * msquared * vars.h[i][j] / vars.chi * Asquared_Re * Avec0_Resquared 
			-      c4 * msquared * vars.h[i][j] / vars.chi * Avec0_Resquared * Avec0_Resquared; 
	 	}

	//////////////////////////Im part /////////////////////////////////
	
	data_t Asquared_Im = 0;	
        const data_t Avec0_Imsquared = vars.Avec0_Im * vars.Avec0_Im ;     	
    	FOR2(i, j)
    	{
        	Asquared_Im += gamma_UU[i][j] * vars.chi * vars.Avec_Im[i] * vars.Avec_Im[j];
    	}

	
	// massive terms 
	rho_potential += 0.5 * msquared * Avec0_Imsquared 
		      + 0.5 * msquared * Asquared_Im; 
	// Self interacting terms 
	rho_potential +=     c4 * msquared * Asquared_Im * Asquared_Im 
		       + 2. * c4 * msquared * Asquared_Im * Avec0_Imsquared 
		       - 3. * c4 * msquared * Avec0_Imsquared * Avec0_Imsquared;

    	FOR1(i)
    	{
				// massive terms 
        	Si_potential[i] += msquared * vars.Avec0_Im * vars.Avec_Im[i]
				// Self interacting terms 
				+ 4. * c4 * msquared * vars.Avec0_Im * vars.Avec_Im[i] * Asquared_Im 
				- 4. * c4 * msquared * vars.Avec0_Im * vars.Avec_Im[i] * Avec0_Imsquared;
    	}

	
    	FOR2(i, j)
    	{
		// massive terms 
        	Sij_potential[i][j] +=
            	msquared * (vars.Avec_Im[i] * vars.Avec_Im[j] +
                        0.5 * vars.h[i][j] / vars.chi * vars.Avec0_Im * vars.Avec0_Im);

            	Sij_potential[i][j] += - 0.5 * vars.h[i][j] / vars.chi *
                                 msquared * Asquared_Im;
        	
    		// Self interacting terms 
        	Sij_potential[i][j] +=
            		  4. * c4 * msquared * vars.Avec_Im[i] * vars.Avec_Im[j] * Asquared_Im
		 	-      c4 * msquared * vars.h[i][j] / vars.chi * Asquared_Im * Asquared_Im
			- 4. * c4 * msquared * vars.Avec_Im[i] * vars.Avec_Im[j] * Avec0_Imsquared
			+ 2. * c4 * msquared * vars.h[i][j] / vars.chi * Asquared_Im * Avec0_Imsquared 
			-      c4 * msquared * vars.h[i][j] / vars.chi * Avec0_Imsquared * Avec0_Imsquared; 
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
        Tensor<2, data_t> DA_Re;
        FOR2(i, j)
        {
            DA_Re[i][j] = d1.Avec_Re[j][i];
            FOR1(k) { DA_Re[i][j] += -chris_phys[k][i][j] * vars.Avec_Re[k]; }
        }

        // DAScalar = D_i A^i
        data_t DA_Re_scalar;
        DA_Re_scalar = 0;
        FOR2(i, j) { DA_Re_scalar += DA_Re[i][j] * gamma_UU[i][j] * vars.chi; }

        // Xsquared = X^/mu X_/mu
        data_t Xsquared_Re;
        Xsquared_Re = -vars.Avec0_Re * vars.Avec0_Re;
        FOR2(i, j) { Xsquared_Re += gamma_UU[i][j] * vars.chi * vars.Avec_Re[j] * vars.Avec_Re[i]; }

        // C = 1 + 4 c4 A^k A_k - 12 c4 Avec0_Re^2
        data_t C_Re = 1.0 - 12.0 * c4 * vars.Avec0_Re * vars.Avec0_Re;
        FOR2(i, j)
        {
            C_Re += 4.0 * c4 * gamma_UU[i][j] * vars.chi * vars.Avec_Re[j] * vars.Avec_Re[i];
        }

        // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
        dVdA_Re = pow(m_params.mass, 2.0) * (1.0 + 4.0 * c4 * Xsquared_Re);

        // dphidt - for now the whole thing is here since it depends mainly
        // on the form of the potential - except the advection term which is in
        // the ProcaField code
        dAvec0dt_Re = 0;
        FOR2(i, j)
        {
            dAvec0dt_Re += -gamma_UU[i][j] * vars.chi * vars.Avec_Re[i] * d1.lapse[j];
        }
        // QUESTION: Should this be lapse * Z / C  or lapse * Z??
        dAvec0dt_Re += -vars.lapse * vars.Zvec_Re / C_Re;
        FOR4(i, j, k, l)
        {
            dAvec0dt_Re += -8.0 * c4 * vars.lapse / C_Re * gamma_UU[i][k] *
                      gamma_UU[j][l] * vars.Avec_Re[i] * vars.Avec_Re[j] * DA_Re[k][l];
        }
        dAvec0dt_Re += vars.lapse / C_Re * (1.0 + 4.0 * c4 * Xsquared_Re) *
                  (vars.K * vars.Avec0_Re - DA_Re_scalar);
        FOR1(i)
        {
            dAvec0dt_Re += 8.0 * c4 * vars.Avec0_Re * vars.lapse / C_Re *
                      (vars.Evec_Re[i] * vars.Avec_Re[i]);
        }
        FOR4(i, j, k, l)
        {
            dAvec0dt_Re += 8.0 * c4 * vars.Avec0_Re * vars.lapse / C_Re *
                      (-K_tensor[i][j] * vars.Avec_Re[k] *
                       vars.Avec_Re[l] * gamma_UU[i][k] * gamma_UU[j][l]);
        }
        FOR2(i, j)
        {
            dAvec0dt_Re += 8.0 * c4 * vars.Avec0_Re * vars.lapse / C_Re *
                      (2.0 * vars.Avec_Re[i] * d1.Avec0_Re[j] * gamma_UU[i][j] * vars.chi);
        }

///////////////////////////Im part////////////////////////////////////

	// Here we are defining often used terms
        // DA[i][j] = D_i A_j
        Tensor<2, data_t> DA_Im;
        FOR2(i, j)
        {
            DA_Im[i][j] = d1.Avec_Im[j][i];
            FOR1(k) { DA_Im[i][j] += -chris_phys[k][i][j] * vars.Avec_Im[k]; }
        }

        // DAScalar = D_i A^i
        data_t DA_Im_scalar;
        DA_Im_scalar = 0;
        FOR2(i, j) { DA_Im_scalar += DA_Im[i][j] * gamma_UU[i][j] * vars.chi; }

        // Xsquared = X^/mu X_/mu
        data_t Xsquared_Im;
        Xsquared_Im = -vars.Avec0_Im * vars.Avec0_Im;
        FOR2(i, j) { Xsquared_Im += gamma_UU[i][j] * vars.chi * vars.Avec_Im[j] * vars.Avec_Im[i]; }

        // C = 1 + 4 c4 A^k A_k - 12 c4 Avec0_Im^2
        data_t C_Im = 1.0 - 12.0 * c4 * vars.Avec0_Im * vars.Avec0_Im;
        FOR2(i, j)
        {
            C_Im += 4.0 * c4 * gamma_UU[i][j] * vars.chi * vars.Avec_Im[j] * vars.Avec_Im[i];
        }

        // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
        dVdA_Im = pow(m_params.mass, 2.0) * (1.0 + 4.0 * c4 * Xsquared_Im);

        // dphidt - for now the whole thing is here since it depends mainly
        // on the form of the potential - except the advection term which is in
        // the ProcaField code
        dAvec0dt_Im = 0;
        FOR2(i, j)
        {
            dAvec0dt_Im += -gamma_UU[i][j] * vars.chi * vars.Avec_Im[i] * d1.lapse[j];
        }
        // QUESTION: Should this be lapse * Z / C  or lapse * Z??
        dAvec0dt_Im += -vars.lapse * vars.Zvec_Im / C_Im;
        FOR4(i, j, k, l)
        {
            dAvec0dt_Im += -8.0 * c4 * vars.lapse / C_Im * gamma_UU[i][k] *
                      gamma_UU[j][l] * vars.Avec_Im[i] * vars.Avec_Im[j] * DA_Im[k][l];
        }
        dAvec0dt_Im += vars.lapse / C_Im * (1.0 + 4.0 * c4 * Xsquared_Im) *
                  (vars.K * vars.Avec0_Im - DA_Im_scalar);
        FOR1(i)
        {
            dAvec0dt_Im += 8.0 * c4 * vars.Avec0_Im * vars.lapse / C_Im *
                      (vars.Evec_Im[i] * vars.Avec_Im[i]);
        }
        FOR4(i, j, k, l)
        {
            dAvec0dt_Im += 8.0 * c4 * vars.Avec0_Im * vars.lapse / C_Im *
                      (-K_tensor[i][j] * vars.Avec_Im[k] *
                       vars.Avec_Im[l] * gamma_UU[i][k] * gamma_UU[j][l]);
        }
        FOR2(i, j)
        {
            dAvec0dt_Im += 8.0 * c4 * vars.Avec0_Im * vars.lapse / C_Im *
                      (2.0 * vars.Avec_Im[i] * d1.Avec0_Im[j] * gamma_UU[i][j] * vars.chi);
        }
    }
};

#endif /* POTENTIAL_HPP_ */
