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
//#include "FixedBGProcaField.hpp"
#include "InitialConditions.hpp"
#include "Potential.hpp"
#include "XSquared.hpp"
#include "ComplexProcaField.hpp"
#include "GammaCalculator.hpp"
#include "NewMatterConstraints.hpp"
#include "ComplexAvecTaggingCriterion.hpp"
#include "Density.hpp"
#include "ComplexProcaFieldConstraints.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"

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

    InitialConditions set_field( m_dx , m_p.initalcondition_data);

    BoxLoops::loop( make_compute_pack(SetValue(0.), set_field), m_state_new, m_state_new, FILL_GHOST_CELLS,disable_simd());

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

}

// Things to do after each timestep
void ProcaFieldLevel::specificPostTimeStep()
{
    CH_TIME("BinaryBHLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    if (m_p.activate_extraction == 1)
    {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(Weyl4(m_p.extraction_params.center, m_dx),
                           m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);

            // Do the extraction on the min extraction level
            if (m_level == min_level)
            {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;
                m_gr_amr.m_interpolator->refresh(fill_ghosts);
                m_gr_amr.fill_multilevel_ghosts(
                    VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                    min_level);
                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);
                my_extraction.execute_query(m_gr_amr.m_interpolator);

            }
        }
    }

}

// Things to do before outputting a plot file
void ProcaFieldLevel::prePlotLevel() {

    fillAllGhosts();
    Potential potential;
    ComplexProcaField<Potential> proca_field(potential, m_p.proca_damping);
    BoxLoops::loop(
        MatterConstraints<ComplexProcaField<Potential>>(
            proca_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    BoxLoops::loop(
        ComplexProcaFieldConstraints<ComplexProcaField<Potential>>(
            proca_field, m_dx, m_p.field_mu, m_p.G_Newton),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    BoxLoops::loop(
        Density<ComplexProcaField<Potential>>(
            proca_field, m_dx),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void ProcaFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                      const double a_time)
{
    // Calculate MatterCCZ4 right hand side with matter_t = ProcaField
    ComplexProcaField<Potential> proca_field(m_p.field_mu, m_p.proca_damping);
    MatterCCZ4<ComplexProcaField<Potential>> my_ccz4_matter(
            proca_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
            m_p.G_Newton);

    BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);

}

void ProcaFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                              const FArrayBox &current_state)
{
    const double radius_bh = 1.75;
    BoxLoops::loop(ComplexAvecTaggingCriterion(m_dx,
		    m_p.threshold_phi,m_p.threshold_K, m_level,m_p.matter_max_level,m_p.extraction_params,m_p.activate_extraction  ),
                   current_state, tagging_criterion, disable_simd());
}
