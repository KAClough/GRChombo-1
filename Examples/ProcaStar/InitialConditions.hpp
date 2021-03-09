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
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "ComplexProcaField.hpp"

//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    const double m_spacing;
    const double m_omega;
    const std::vector<double> m_a0;
    const std::vector<double> m_da0dr;
    const std::vector<double> m_a1;
    const std::vector<double> m_m;
    const std::vector<double> m_sig;
    const std::array<double, CH_SPACEDIM> m_center;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t>
    using Vars = ComplexProcaField::Vars<data_t>;

  public:
    struct params_t {
        std::vector<double> a0;
        std::vector<double> da0dr;
        std::vector<double> a1;
        std::vector<double> m;
        std::vector<double> sig;
        double omega;
        double spacing;
    };


    //! The constructor for the class
    InitialConditions(
                      const std::array<double, CH_SPACEDIM> a_center,
                      const double a_dx,
                      params_t initalcondition_data
                      )
        : m_dx(a_dx), m_center(a_center), m_a0(initalcondition_data.a0),
        m_da0dr(initalcondition_data.da0dr),m_a1(initalcondition_data.a1),
        m_m(initalcondition_data.m),m_sig(initalcondition_data.sig),
        m_spacing(initalcondition_data.spacing),m_omega(initalcondition_data.omega)
    {
    }


     double linear_interpolation(const  std::vector<double> vector, const double rr) const {
        const int indxL = static_cast<int>(floor(rr / m_spacing));
        const int indxH = static_cast<int>(ceil(rr / m_spacing));
        const int ind_max = vector.size();

        if ( indxH < ind_max){
            const double interpolation_value =
            vector[indxL] +
            (rr / m_spacing - indxL) * (vector[indxH] - vector[indxL]);

            return interpolation_value;
        }
        else{
            return vector[ind_max-1];
        }
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        Vars<data_t> vars;
        VarsTools::assign(vars, 0.);

        Tensor<2,double> g; // Metix Index low low
        Tensor<2,double> g_conf; // Metix Index low low
        Tensor<2,double> g_spher; // Metix Index low low
        Tensor<2,double> jacobian; // Metix Index low low
        Tensor<1,double> Avec_spher_Re;
        Tensor<1,double> Avec_spher_Im;
        Tensor<1,double> Avec_Re;
        Tensor<1,double> Avec_Im;
        Tensor<1,double> Evec_spher_Re;
        Tensor<1,double> Evec_spher_Im;
        Tensor<1,double> Evec_Re;
        Tensor<1,double> Evec_Im;
        Tensor<1,double> Evec_Re_U;
        Tensor<1,double> Evec_Im_U;
        FOR2(i,j){ g_spher[i][j] = 0; g[i][j] = 0;}
        FOR1(i){
                 Avec_spher_Re[i] = 0;
                 Avec_spher_Im[i] = 0;
                 Avec_Re[i] = 0;
                 Avec_Im[i] = 0;
                 Evec_spher_Re[i] = 0;
                 Evec_spher_Im[i] = 0;
                 Evec_Re[i] = 0;
                 Evec_Im[i] = 0;
                 Evec_Re_U[i] = 0;
                 Evec_Im_U[i] = 0;
                }
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double t = 0;

        const double rr = coords.get_radius();
        const double rr2 = rr*rr;
        double rho2 = pow(x, 2) + pow(y, 2);

        double rho = sqrt(rho2);
        // sinus(theta):
        double sintheta = rho / rr;
        double costheta = z / rr;
        // cos(phi)
        double cosphi = x / rho;
        // sin(phi)
        double sinphi = y / rho;

        const double norm = 1.0/sqrt(16.0*M_PI);
        const double a0 = norm*linear_interpolation(m_a0,rr);
        const double da0dr = norm*linear_interpolation(m_da0dr,rr);
        const double a1 = norm*linear_interpolation(m_a1,rr);
        const double m = linear_interpolation(m_m,rr);
        const double sig = linear_interpolation(m_sig,rr);

        jacobian[0][0] = x/rr ;
        jacobian[1][0] = cosphi*z/rr2 ;
        jacobian[2][0] = -y/rho2 ;
        jacobian[0][1] = y/rr ;
        jacobian[1][1] = sinphi*z/rr2 ;
        jacobian[2][1] = x/rho2 ;
        jacobian[0][2] = z/rr ;
        jacobian[1][2] = -rho/rr2 ;
        jacobian[2][2] = 0.0 ;

        const double lapse = sig * sqrt(1.0 - 2.0* m / rr);
        // Define r r component
        g_spher[0][0] = 1.0/(1.0 - 2.0 * m / rr );
        // Define theta theta component
        g_spher[1][1] = rr2 ;
        // Define phi phi component
        g_spher[2][2] = rr2 * pow(sintheta, 2);

        // set the field variable to approx profile
        data_t phi_Re = - 1.0/lapse * a0 * cos(-m_omega*t);
        data_t phi_Im = - 1.0/lapse * a0 * sin(-m_omega*t);
        // r Component
        Avec_spher_Re[0] =  a1 * sin(m_omega * t);
        Avec_spher_Im[0] =  a1 * cos(m_omega * t);
        Evec_spher_Re[0] =  1.0/(sqrt(1.0 - 2.0 * m / rr )*sig) * (-m_omega*a1 + da0dr) * cos(-m_omega * t);
        Evec_spher_Im[0] =  1.0/(sqrt(1.0 - 2.0 * m / rr )*sig) * (-m_omega*a1 + da0dr) * sin(-m_omega * t);

        FOR2(i, j)
        {
            Avec_Re[i] += Avec_spher_Re[j] * jacobian[j][i];
            Avec_Im[i] += Avec_spher_Im[j] * jacobian[j][i];
            Evec_Re[i] += Evec_spher_Re[j] * jacobian[j][i];
            Evec_Im[i] += Evec_spher_Im[j] * jacobian[j][i];
            Evec_Re_Ref[i] += Evec_spher_Re_Ref[j] * jacobian_inverse[j][i];
            Evec_Im_Ref[i] += Evec_spher_Im_Ref[j] * jacobian_inverse[j][i];
            FOR2(k, l)
                    {
                        g[i][j] += g_spher[k][l] * jacobian[k][i] * jacobian[l][j];
                    }
        }

        const auto g_UU = TensorAlgebra::compute_inverse(g);

        FOR2(i, j)
        {
            Evec_Re_U[i] += Evec_Re[j] * g_UU[j][i];
            Evec_Im_U[i] += Evec_Im[j] * g_UU[j][i];
        }
         data_t deth = TensorAlgebra::compute_determinant<data_t>(g);

         data_t chi = pow(deth,-1.0/3.0);

         FOR2(i,j)
         {
           g_conf[i][j] = g[i][j]*chi;
         }


        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(lapse, c_lapse);

        current_cell.store_vars(g[0][0], c_h11);
        current_cell.store_vars(g[0][1], c_h12);
        current_cell.store_vars(g[0][2], c_h13);
        current_cell.store_vars(g[1][1], c_h22);
        current_cell.store_vars(g[1][2], c_h23);
        current_cell.store_vars(g[2][2], c_h33);

        current_cell.store_vars(Avec_Re[0], c_Avec1_Re);
        current_cell.store_vars(Avec_Re[1], c_Avec2_Re);
        current_cell.store_vars(Avec_Re[2], c_Avec3_Re);

        current_cell.store_vars(Avec_Im[0], c_Avec1_Im);
        current_cell.store_vars(Avec_Im[1], c_Avec2_Im);
        current_cell.store_vars(Avec_Im[2], c_Avec3_Im);

        current_cell.store_vars(Evec_Re_U[0], c_Evec1_Re);
        current_cell.store_vars(Evec_Re_U[1], c_Evec2_Re);
        current_cell.store_vars(Evec_Re_U[2], c_Evec3_Re);

        current_cell.store_vars(Evec_Im_U[0], c_Evec1_Im);
        current_cell.store_vars(Evec_Im_U[1], c_Evec2_Im);
        current_cell.store_vars(Evec_Im_U[2], c_Evec3_Im);

        current_cell.store_vars(phi_Re, c_Avec0_Re);
        current_cell.store_vars(phi_Im, c_Avec0_Im);

    }
};

#endif /* INITIALCONDITIONS_HPP_ */
