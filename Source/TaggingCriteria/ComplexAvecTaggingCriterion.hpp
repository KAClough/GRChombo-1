/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXAVECTAGGINGCRITERION_HPP_
#define COMPLEXAVECTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class ComplexAvecTaggingCriterion
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_phi;
    const double m_threshold_K;
    const int m_max_matter_level;
    const int m_level;
    const SphericalExtraction::params_t m_params;
    const bool m_activate_extraction;
	
  public:
    ComplexAvecTaggingCriterion(double dx, double threshold_phi,
			       double threshold_K, int a_level, 
			       int max_matter_level, 
			       const SphericalExtraction::params_t a_params,
                               const bool activate_extraction = false)
        : m_dx(dx), m_deriv(dx), m_threshold_phi(threshold_phi),
          m_threshold_K(threshold_K), m_level(a_level), 
	  m_max_matter_level(max_matter_level), m_params(a_params),
          m_activate_extraction(activate_extraction)
	 {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Tensor<1, data_t> d1_phi_re;
        FOR1(idir) m_deriv.diff1(d1_phi_re, current_cell, idir, c_Avec0_Re);

        Tensor<1, data_t> d1_phi_im;
        FOR1(idir) m_deriv.diff1(d1_phi_im, current_cell, idir, c_Avec0_Im);

        Tensor<1, data_t> d1_K;
        FOR1(idir) m_deriv.diff1(d1_K, current_cell, idir, c_K);

        data_t mod_d1_phi = 0;
        data_t mod_d1_K = 0;
        FOR1(idir)
        {
	    if(m_level < m_max_matter_level){
            mod_d1_phi += d1_phi_re[idir] * d1_phi_re[idir] +
                          d1_phi_im[idir] * d1_phi_im[idir];
	    } 
            mod_d1_K += d1_K[idir] * d1_K[idir];
        }

        data_t criterion = m_dx * (1.0 / m_threshold_phi * sqrt(mod_d1_phi) +
                                   1.0 / m_threshold_K * sqrt(mod_d1_K));

        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
 
                }
            }
	}

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* COMPLEXAVECTAGGINGCRITERION_HPP_ */
