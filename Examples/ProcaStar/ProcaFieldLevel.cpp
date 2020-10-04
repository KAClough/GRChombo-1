/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ProcaFieldLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "NanCheck.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"

// For tag cells
#include "FixedGridsTaggingCriterionBH.hpp"
#include "MatterCCZ4.hpp"

// Problem specific includes
#include "ExcisionProcaDiagnostics.hpp"
#include "ExcisionProcaEvolution.hpp"
#include "FixedBGDensityAndAngularMom.hpp"
#include "FixedBGEnergyAndAngularMomFlux.hpp"
#include "FixedBGEvolution.hpp"
//#include "FixedBGProcaField.hpp"
#include "FluxExtraction.hpp"
#include "InitialConditions.hpp"
#include "KerrSchildFixedBG.hpp"
#include "Potential.hpp"
#include "XSquared.hpp"
#include "ComplexProcaField.hpp"
#include "GammaCalculator.hpp"
#include "NewMatterConstraints.hpp"
//#include "ProcaConstraint.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ProcaFieldLevel::specificAdvance()
{
    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(m_dx), m_state_new, m_state_new, SKIP_GHOST_CELLS,
                       disable_simd());
}

// Initial data for field and metric variables
void ProcaFieldLevel::initialData()
{
    CH_TIME("ProcaFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ProcaFieldLevel::initialData " << m_level << endl;

    // Set the ICs

    InitialConditions set_field( m_p.center,
                               m_dx , m_p.a0, m_p.a1, m_p.m,
                              m_p.sig, m_p.spacing, m_p.omega);

    BoxLoops::loop(set_field, m_state_new, m_state_new, FILL_GHOST_CELLS,disable_simd());

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);
    // BoxLoops::loop(kerr_bh, m_state_new, m_state_diagnostics,
    //                SKIP_GHOST_CELLS);

    // now the gauss constraint
    /*
        fillAllGhosts();
        ProcaConstraint enforce_constraint(m_p.center, m_p.bg_params,
                                           m_p.potential_params.mass, m_dx);
        BoxLoops::loop(enforce_constraint, m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS);
    */
    // make excision data zero
}

// Things to do after each timestep
void ProcaFieldLevel::specificPostTimeStep()
{}

// Things to do before outputting a plot file
void ProcaFieldLevel::prePlotLevel() {

    fillAllGhosts();
    ComplexProcaField proca_field(m_p.field_mu, m_p.proca_damping);
    BoxLoops::loop(
        MatterConstraints<ComplexProcaField>(
            proca_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

}

// Things to do in RHS update, at each RK4 step
void ProcaFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                      const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ProcaField
    ComplexProcaField proca_field(m_p.field_mu, m_p.proca_damping);
    MatterCCZ4<ComplexProcaField> my_ccz4_matter(
            proca_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
            m_p.G_Newton);

    BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);

}

void ProcaFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                              const FArrayBox &current_state)
{
    const double radius_bh = 1.75;
    BoxLoops::loop(FixedGridsTaggingCriterionBH(m_dx, m_level, m_p.max_level, m_p.L, m_p.center, radius_bh),
                   current_state, tagging_criterion, disable_simd());
}
