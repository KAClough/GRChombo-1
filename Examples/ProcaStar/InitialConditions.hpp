/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "KerrSchildFixedBG.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    const double m_amplitude;
    const double m_spacing;
    const double m_omega;
    const std::vector<double> m_a0;
    const std::vector<double> m_a1;
    const std::vector<double> m_m;
    const std::vector<double> m_sig;
    const std::array<double, CH_SPACEDIM> m_center;
    const KerrSchildFixedBG::params_t m_bg_params;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  public:
    //! The constructor for the class
    InitialConditions(const double a_amplitude,
                      const std::array<double, CH_SPACEDIM> a_center,
                      const KerrSchildFixedBG::params_t a_bg_params,
                      const double a_dx,
                      std::vector<double> a0,
                      std::vector<double> a1,
                      std::vector<double> m,
                      std::vector<double> sig,
                      double spacing,
                      double omega
                      )
        : m_dx(a_dx), m_amplitude(a_amplitude), m_center(a_center),
          m_bg_params(a_bg_params), m_a0(a0), m_a1(a1), m_m(m), m_sig(sig),
          m_spacing(spacing), m_omega(omega)
    {
    }


     double linear_interpolation(const  std::vector<double> vector, const double rr) const {
        int indxL = static_cast<int>(floor(rr / m_spacing));
        int indxH = static_cast<int>(ceil(rr / m_spacing));

        if ( indxH < vector.size()){
            const double interpolation_value =
            vector[indxL] +
            (rr / m_spacing - indxL) * (vector[indxH] - vector[indxL]);

            return interpolation_value;
        }
        else{

            return 0;
        }
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        Tensor<2,double> g; // Metix Index low low
        Tensor<2,double> g_spher; // Metix Index low low
        Tensor<2,double> jacobian; // Metix Index low low
        Tensor<1,double> Avec_spher;
        Tensor<1,double> Avec;
        FOR2(i,j){ g_spher[i][j] = 0; g[i][j] = 0;}
        FOR1(i){ Avec_spher[i] = 0;}
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double t = 0;

        const double rr = coords.get_radius();
        const double rr2 = rr*rr;
        double rho2 = pow(x, 2) + pow(y, 2);

        double rho = sqrt(rho2);
        // sinus(theta)
        double sintheta = rho / rr;
        double costheta = z / rr;
        // cos(phi)
        double cosphi = x / rho;
        // sin(phi)
        double sinphi = y / rho;


        const double a0 = linear_interpolation(m_a0,rr);
        const double a1 = linear_interpolation(m_a1,rr);
        const double m = linear_interpolation(m_m,rr);
        const double sig = linear_interpolation(m_sig,rr);

        jacobian[0][0] = x / rr;
        jacobian[1][0] = cosphi * z / rr2;
        jacobian[2][0] = -y / rho2;
        jacobian[0][1] = y / rr;
        jacobian[1][1] = sinphi * z / rr2;
        jacobian[2][1] = x / rho2;
        jacobian[0][2] = z / rr;
        jacobian[1][2] = -rho / rr2;
        jacobian[2][2] = 0.0;
        jacobian[3][3] = 1.0;

        const double lapse = sig * sqrt(1.0 - 2.0* m / rr);
        // Define r r component
        g_spher[0][0] = 1.0/(1.0 - 2.0 * m / rr );
        // Define theta theta component
        g_spher[1][1] = rr2 ;
        // Define phi phi component
        g_spher[2][2] = rr2 * pow(sintheta, 2);

        // set the field variable to approx profile
        data_t phi = a0*cos(m_omega*t);
        // r Component
        Avec_spher[0] = - a1 * sin(m_omega * t);

        FOR2(i, j)
        {
            Avec[i] += Avec_spher[i] * jacobian[i][j];
            FOR2(k, l)
                    {
                        g[i][j] += g_spher[k][l] * jacobian[k][i] * jacobian[l][j];
                    }
        }



        current_cell.store_vars(phi, c_phi);
    }
};

#endif /* INITIALCONDITIONS_HPP_ */
