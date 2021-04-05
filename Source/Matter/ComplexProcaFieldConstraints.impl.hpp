/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(COMPLEXPROCAFIELDCONSTRAINTS_HPP_)
#error "This file should only be included through ComplexProcaFieldConstraints.hpp"
#endif

#ifndef COMPLEXPROCAFIELDCONSTRAINTS_IMPL_HPP_
#define COMPLEXPROCAFIELDCONSTRAINTS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
ComplexProcaFieldConstraints<matter_t>::ComplexProcaFieldConstraints(const matter_t a_matter,
                                               double dx, double vector_mass, double G_Newton)
    : my_matter(a_matter), m_G_Newton(G_Newton), m_deriv(dx), m_vector_mass(vector_mass)
{
}

template <class matter_t>
template <class data_t>
void ComplexProcaFieldConstraints<matter_t>::compute(Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto d1 = m_deriv.template diff1<BSSNMatterVars>(current_cell);
    const auto d2 = m_deriv.template diff2<BSSNMatterVars>(current_cell);
    const auto vars = current_cell.template load_vars<BSSNMatterVars>();


    // Inverse metric and Christoffel symbol
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);
    const Tensor<3, data_t> chris_phys = TensorAlgebra::compute_phys_chris(
        d1.chi, vars.chi, vars.h, h_UU, chris.ULL);
    // Energy Momentum Tensor
    const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    // Hamiltonian constraint
    data_t Proca = pow(m_vector_mass, 2.0) * vars.Avec0_Re;

    FOR1(i)
    {
        Proca +=  d1.Evec_Re[i][i];
        FOR1(j)
        {
            Proca +=
                chris_phys[i][i][j] * vars.Evec_Re[j];
        }
    }
    // Write the constraints into the output FArrayBox
     current_cell.store_vars(Proca, c_Proca);
}

#endif /* COMPLEXPROCAFIELDCONSTRAINTS_IMPL_HPP_ */
