/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef VIOLATION_HPP_
#define VIOLATION_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FixedBGProcaField.hpp"
#include "IsotropicKerrFixedBG.hpp"
#include "Potential.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class violation
{
  protected:
    const Potential::params_t m_potential_params;
    const IsotropicKerrFixedBG::params_t m_bg_params;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables


    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t>
    using Vars = FixedBGProcaField<Potential>::template Vars<data_t>;

  public:
    //! The constructor for the class
    violation(const Potential::params_t a_potential_params,
             const IsotropicKerrFixedBG::params_t a_bg_params,
             const std::array<double, CH_SPACEDIM> a_center, const double a_dx)
        : m_dx(a_dx), m_center(a_center),
          m_potential_params(a_potential_params), m_bg_params(a_bg_params),
	  m_deriv(a_dx)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        IsotropicKerrFixedBG kerr_bh(m_bg_params, m_dx);
        MetricVars<data_t> metric_vars;
        kerr_bh.compute_metric_background(metric_vars, current_cell);
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        // The matter vars
        const auto vars = current_cell.template load_vars<Vars>();

	// calculate full spatial christoffel symbols
    	using namespace TensorAlgebra;
    	const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        const auto d1 = m_deriv.template diff1<Vars>(current_cell);

        // Xsquared = X^/mu X_/mu
        data_t violation;

	Potential massive_photon(m_potential_params, m_bg_params );
	
	//calculate r_plus values
	double r_plus = m_bg_params.mass + pow( m_bg_params.mass * m_bg_params.mass - 
			m_bg_params.spin * m_bg_params.spin, 0.5);

        violation = metric_vars.lapse * (pow( massive_photon.m_field(coords.x,coords.y,coords.z,r_plus),
			       	2.0) * vars.phi);

	FOR1(i)
	    {
		violation += metric_vars.lapse * d1.Evec[i][i];
		FOR1(j)
		{
		    violation += metric_vars.lapse * chris_phys.ULL[i][i][j] * vars.Evec[j];
		}
	    }
        // Store the initial values of the variables
        current_cell.store_vars(violation, c_violation);
    }
};

#endif /* VIOLATION_HPP_ */
