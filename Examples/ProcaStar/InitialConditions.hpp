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
#include "Potential.hpp"

//! Class which creates the initial conditions
class InitialConditions
{
  protected:
    const double m_dx;
    const double m_spacing;
    const double m_omega;
    const double* m_a0;
    const double* m_da0dr;
    const double* m_a1;
    const double* m_m;
    const double* m_sig;
    const int m_size;
    const std::array<double, CH_SPACEDIM> m_center1;
    const std::array<double, CH_SPACEDIM> m_center2;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    // The evolution vars
    template <class data_t>
    using Vars = ComplexProcaField<Potential>::Vars<data_t>;

  public:
    struct params_t {
        double* a0;
        double* da0dr;
        double* a1;
        double* m;
        double* sig;
	int size;
        double omega;
        double spacing;
    	std::array<double, CH_SPACEDIM> center1;
    	std::array<double, CH_SPACEDIM> center2;
    };
    
    struct proca_star_struct {
	double lapse;
	double phi_Re;
	double phi_Im;
        Tensor<2,double> g; // Metric Index low low
        Tensor<1,double> Avec_Re;
        Tensor<1,double> Avec_Im;
        Tensor<1,double> Evec_Re_U;
        Tensor<1,double> Evec_Im_U;
    };

    //! The constructor for the class
    InitialConditions(const double a_dx,
                      const params_t initalcondition_data
                      )
        : m_dx(a_dx), m_center1(initalcondition_data.center1),  m_center2(initalcondition_data.center2),
	m_a0(initalcondition_data.a0), m_da0dr(initalcondition_data.da0dr), m_a1(initalcondition_data.a1),
        m_m(initalcondition_data.m),m_sig(initalcondition_data.sig),
        m_spacing(initalcondition_data.spacing),m_omega(initalcondition_data.omega), m_size(initalcondition_data.size)
    {
    }


     double linear_interpolation(const double* vector, const double rr) const {
        const int indxL = static_cast<int>(floor(rr / m_spacing));
        const int indxH = static_cast<int>(ceil(rr / m_spacing));
        const int ind_max = m_size;

        if ( indxH < ind_max){
            const double interpolation_value =
            *(vector+indxL) +
            (rr / m_spacing - indxL) * (*(vector+indxH) - *(vector+indxL));

            return interpolation_value;
        }
        else{
            return *(vector+ind_max-1);
        }
    }
    

    template <class data_t> void get_proca_star_values( const double x, const double y, const double z, const double t,  proca_star_struct &proca_star_values ) const {

        Tensor<2,double> g; // Metric Index low low
        Tensor<2,double> g_spher; // Metric Index low low
        Tensor<2,double> jacobian; // Metric Index low low
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

        const double rr = sqrt(x*x+y*y+z*z);
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

        const double norm = 1.0/sqrt(16*M_PI);
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
        double phi_Re = - 1.0/lapse * a0 * cos(m_omega*t);
        double phi_Im = - 1.0/lapse * a0 * sin(m_omega*t);
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
            FOR2(k, l)
                    {
                        g[i][j] += g_spher[k][l] * jacobian[k][i] * jacobian[l][j];
                    }
        }

         const auto g_UU = TensorAlgebra::compute_inverse(g);

         // Raising the electric vector

         FOR2(i, j)
         {
            Evec_Re_U[i] += Evec_Re[j] * g_UU[j][i];
            Evec_Im_U[i] += Evec_Im[j] * g_UU[j][i];
         }
	proca_star_values.phi_Re = phi_Re;
	proca_star_values.phi_Im = phi_Im;
	proca_star_values.lapse = lapse;
	proca_star_values.g = g; // Metric Index low low
        proca_star_values.Avec_Re = Avec_Re ;
        proca_star_values.Avec_Im = Avec_Im;
        proca_star_values.Evec_Re_U = Evec_Re_U;
        proca_star_values.Evec_Im_U = Evec_Im_U;
    }
 

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Vars<data_t> vars;
        VarsTools::assign(vars, 0.);
        Tensor<2,double> g_conf; // Metric Index low low

	// Get first star data 
        Coordinates<data_t> coords1(current_cell, m_dx, m_center1);
        const double x1 = coords1.x;
        const double y1 = coords1.y;
        const double z1 = coords1.z;
        const double t1 = 0;

	proca_star_struct star1;
	get_proca_star_values<data_t>(x1,y1,z1,t1,star1);

	// Get second star data 
        Coordinates<data_t> coords2(current_cell, m_dx, m_center2);
        const double x2 = coords2.x;
        const double y2 = coords2.y;
        const double z2 = coords2.z;
        const double t2 = 0;

	proca_star_struct star2;
	get_proca_star_values<data_t>(x2,y2,z2,t2,star2);

	// Get corrections 
        const double x_corr = (m_center2[0] - m_center1[0]);
        const double y_corr = (m_center2[1] - m_center1[1]);
        const double z_corr = (m_center2[2] - m_center1[2]);
        const double t_corr = 0;


	proca_star_struct corr;
	get_proca_star_values<data_t>(x_corr,y_corr,z_corr,t_corr,corr);

	// Superpose two stars 	
        Tensor<1,double> Avec_Re;
        Tensor<1,double> Avec_Im;
        Tensor<1,double> Evec_Re_U;
        Tensor<1,double> Evec_Im_U;
        Tensor<2,double> g; // Metric Index low low
	Tensor<2,double> h;

	// Secret sauce fix 
	FOR2(i,j) h[i][j] = corr.g[i][j]; 

	/*FOR2(i,j){
		h[i][j] = 0 ; 
	}
	FOR1(i) h[i][i] = 1;
	*/
	double lapse = star1.lapse + star2.lapse - 1.0 ;
	double phi_Re = star1.phi_Re + star2.phi_Re;
	double phi_Im = star1.phi_Im + star2.phi_Im;
	FOR1(i){
        	Avec_Re[i] = star1.Avec_Re[i] + star2.Avec_Re[i];
        	Avec_Im[i] = star1.Avec_Im[i] + star2.Avec_Im[i];
        	Evec_Re_U[i] = star1.Evec_Re_U[i] + star2.Evec_Re_U[i];
        	Evec_Im_U[i] = star1.Evec_Im_U[i] + star2.Evec_Im_U[i];
		FOR1(j){
        		g[i][j] = star1.g[i][j]  + star2.g[i][j] - h[i][j]; // Metric Index low low
		}
	 }

         data_t deth = TensorAlgebra::compute_determinant<data_t>(g);

         data_t chi = pow(deth,-1.0/3.0);

         FOR2(i,j)
         {
           g_conf[i][j] = g[i][j]*chi;
         }


        current_cell.store_vars(chi, c_chi);
        current_cell.store_vars(lapse, c_lapse);

        current_cell.store_vars(g_conf[0][0], c_h11);
        current_cell.store_vars(g_conf[0][1], c_h12);
        current_cell.store_vars(g_conf[0][2], c_h13);
        current_cell.store_vars(g_conf[1][1], c_h22);
        current_cell.store_vars(g_conf[1][2], c_h23);
        current_cell.store_vars(g_conf[2][2], c_h33);

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
