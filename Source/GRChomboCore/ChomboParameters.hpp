/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CHOMBOPARAMETERS_HPP_
#define CHOMBOPARAMETERS_HPP_

// General includes
#include "BoundaryConditions.hpp"
#include "GRParmParse.hpp"
#include "UserVariables.hpp"
#include "VariableType.hpp"
#include <algorithm>

class ChomboParameters
{
  public:
    ChomboParameters(GRParmParse &pp) { read_params(pp); }

    void read_params(GRParmParse &pp)
    {
        pp.load("verbosity", verbosity, 0);
        // Grid setup
        pp.load("L", L, 1.0);
        pp.load("regrid_threshold", regrid_threshold, 0.5);
        pp.load("num_ghosts", num_ghosts, 3);
        pp.load("tag_buffer_size", tag_buffer_size, 3);
        pp.load("dt_multiplier", dt_multiplier, 0.25);
        pp.load("fill_ratio", fill_ratio, 0.7);

        // Setup the grid size
        ivN = IntVect::Unit;
        int max_N = 0;
        for (int dir = 0; dir < CH_SPACEDIM; ++dir)
        {
            char dir_str[20];
            sprintf(dir_str, "N%d", dir + 1);
            int N;
            pp.load(dir_str, N);
            N_vect.push_back(N);
            ivN[dir] = N - 1;
            max_N = max(N, max_N);
        }
        coarsest_dx = L / max_N;

        pp.load("center", center,
                {0.5 * N_vect[0] * coarsest_dx, 0.5 * N_vect[1] * coarsest_dx,
                 0.5 * N_vect[2] * coarsest_dx}); // default to center

        // Periodicity and boundaries
        pp.load("isPeriodic", isPeriodic, {true, true, true});
        read_boundary_params(pp);
        pout()
            << "After applying boundary conditions, center has been set to: ";
        FOR1(idir) { pout() << center[idir] << " "; }
        pout() << endl;

        // Misc
        pp.load("ignore_checkpoint_name_mismatch",
                ignore_checkpoint_name_mismatch, false);

        pp.load("max_level", max_level, 0);
        // the reference ratio is hard coded to 2 on all levels
        // in principle it can be set to other values, but this is
        // not recommended since we do not test GRChombo with other
        // refinement ratios - use other values at your own risk
        ref_ratios.resize(max_level + 1);
        ref_ratios.assign(2);
        // read in frequency of regrid on each levels, needs
        // max_level + 1 entries (although never regrids on max_level+1)
        pp.getarr("regrid_interval", regrid_interval, 0, max_level + 1);

        // time stepping outputs and regrid data
        pp.load("checkpoint_interval", checkpoint_interval, 1);
        pp.load("chk_prefix", checkpoint_prefix);
        pp.load("plot_interval", plot_interval, 0);
        pp.load("plot_prefix", plot_prefix);
        pp.load("stop_time", stop_time, 1.0);
        pp.load("max_steps", max_steps, 1000000);
        pp.load("write_plot_ghosts", write_plot_ghosts, false);

        // load vars to write to plot files
        load_vars_to_vector(pp, "plot_vars", "num_plot_vars", plot_vars,
                            num_plot_vars);

        // alias the weird chombo names to something more descriptive
        // for these box params, and default to some reasonable values
        if (pp.contains("max_grid_size"))
        {
            pp.load("max_grid_size", max_grid_size);
        }
        else
        {
            pp.load("max_box_size", max_grid_size, 64);
        }
        if (pp.contains("block_factor"))
        {
            pp.load("block_factor", block_factor);
        }
        else
        {
            pp.load("min_box_size", block_factor, 8);
        }
    }

    // read in of the boundary conditions - this is quite long
    // so kept separate for readability
    void read_boundary_params(GRParmParse &pp)
    {
        // default to static BCs if unspecified in params
        int bc = BoundaryConditions::STATIC_BC;
        pp.load("hi_boundary", boundary_params.hi_boundary, {bc, bc, bc});
        pp.load("lo_boundary", boundary_params.lo_boundary, {bc, bc, bc});
        // set defaults, then override them where appropriate
        boundary_params.vars_parity.fill(BoundaryConditions::EVEN);
        boundary_params.vars_asymptotic_values.fill(0.0);
        boundary_params.is_periodic.fill(true);
        boundary_params.extrapolation_order = 0;
        nonperiodic_boundaries_exist = false;
        boundary_solution_enforced = false;
        boundary_params.mixed_bc_extrapolating_vars.clear();
        boundary_params.mixed_bc_sommerfeld_vars.clear();
        FOR1(idir)
        {
            if (isPeriodic[idir] == false)
            {
                nonperiodic_boundaries_exist = true;
                boundary_params.is_periodic[idir] = false;

                // read in relevent params - note that no defaults are set so as
                // to force the user to specify them where the relevant BCs are
                // selected
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::REFLECTIVE_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::REFLECTIVE_BC))
                {
                    boundary_solution_enforced = true;
                    pp.load("vars_parity", boundary_params.vars_parity);

                    if ((boundary_params.hi_boundary[idir] ==
                         BoundaryConditions::REFLECTIVE_BC) &&
                        (boundary_params.lo_boundary[idir] !=
                         BoundaryConditions::REFLECTIVE_BC))
                    {
                        center[idir] = N_vect[idir] * coarsest_dx;
                    }

                    if ((boundary_params.lo_boundary[idir] ==
                         BoundaryConditions::REFLECTIVE_BC) &&
                        (boundary_params.hi_boundary[idir] !=
                         BoundaryConditions::REFLECTIVE_BC))
                    {
                        center[idir] = 0;
                    }
                }
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::EXTRAPOLATING_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::EXTRAPOLATING_BC) ||
                    (boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::MIXED_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::MIXED_BC))
                {
                    boundary_solution_enforced = true;
                    pp.load("extrapolation_order",
                            boundary_params.extrapolation_order, 1);
                }
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::SOMMERFELD_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::SOMMERFELD_BC) ||
                    (boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::MIXED_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::MIXED_BC))
                {
                    int num_values = -1;
                    std::vector<std::pair<int, VariableType>>
                        nonzero_asymptotic_vars;
                    load_vars_to_vector(pp, "nonzero_asymptotic_vars",
                                        "num_nonzero_asymptotic_vars",
                                        nonzero_asymptotic_vars, num_values);
                    const double default_value = 0.0;
                    load_values_to_array(pp, "nonzero_asymptotic_values",
                                         nonzero_asymptotic_vars,
                                         boundary_params.vars_asymptotic_values,
                                         default_value);

                    // for backwards compatibility, but above method should
                    // be preferred in future
                    if (num_values == -1)
                    {
                        pp.load("vars_asymptotic_values",
                                boundary_params.vars_asymptotic_values);
                    }
                }
                if ((boundary_params.hi_boundary[idir] ==
                     BoundaryConditions::MIXED_BC) ||
                    (boundary_params.lo_boundary[idir] ==
                     BoundaryConditions::MIXED_BC))
                {
                    // only read in these params once, even if specified in
                    // several directions. They must be the same for all.
                    if (boundary_params.mixed_bc_extrapolating_vars.size() > 0)
                        break;

                    // if not yet done, read in the mixed conditions
                    int num_extrapolating_vars = 0;
                    std::vector<std::pair<int, VariableType>>
                        extrapolating_vars;
                    load_vars_to_vector(
                        pp, "extrapolating_vars", "num_extrapolating_vars",
                        extrapolating_vars, num_extrapolating_vars);
                    for (int icomp = 0; icomp < NUM_VARS; icomp++)
                    {
                        bool is_extrapolating = false;
                        // if the variable is not in extrapolating vars, it is
                        // assumed to be sommerfeld by default
                        for (int icomp2 = 0; icomp2 < extrapolating_vars.size();
                             icomp2++)
                        {
                            if (icomp == extrapolating_vars[icomp2].first)
                            {
                                // should be an evolution variable
                                CH_assert(extrapolating_vars[icomp2].second ==
                                          VariableType::evolution);
                                boundary_params.mixed_bc_extrapolating_vars
                                    .push_back(icomp);
                                is_extrapolating = true;
                            }
                        }
                        if (!is_extrapolating)
                        {
                            boundary_params.mixed_bc_sommerfeld_vars
                                .push_back(icomp);
                        }
                    }
                }
            }
        }

        // write out boundary conditions where non periodic - useful for debug
        if (nonperiodic_boundaries_exist)
        {
            BoundaryConditions::write_boundary_conditions(boundary_params);
        }
    }

    // function to create a vector of enums of vars by reading in their
    // names as strings from the params file and converting it to the enums
    void load_vars_to_vector(
        GRParmParse &pp, const char *a_vars_vector_string,
        const char *a_vector_size_string,
        std::vector<std::pair<int, VariableType>> &a_vars_vector,
        int &a_vars_vector_size)
    {
        int num_values;
        pp.load(a_vector_size_string, num_values, -1);
        // only set a_vars_vector and a_var_vector_size if a_vector_size_string
        // found
        if (num_values >= 0)
        {
            std::vector<std::string> var_names(num_values, "");
            pp.load(a_vars_vector_string, var_names, num_values, var_names);
            for (std::string var_name : var_names)
            {
                // first assume plot_var is a normal evolution var
                int var = UserVariables::variable_name_to_enum(var_name);
                VariableType var_type = VariableType::evolution;
                if (var < 0)
                {
                    // if not an evolution var check if it's a diagnostic var
                    var = DiagnosticVariables::variable_name_to_enum(var_name);
                    if (var < 0)
                    {
                        // it's neither :(
                        pout() << "Variable with name " << var_name
                               << " not found." << endl;
                    }
                    else
                    {
                        var_type = VariableType::diagnostic;
                    }
                }
                if (var >= 0)
                {
                    a_vars_vector.emplace_back(var, var_type);
                }
            }
            // overwrites read in value if entries have been ignored
            a_vars_vector_size = a_vars_vector.size();
        }
    }

    // where one has read in a subset of variables with some feature
    // this reads in a set of associated values and assigns it into a full
    // array of all NUM_VARS vars (setting other values to a default value)
    template <class T>
    void load_values_to_array(
        GRParmParse &pp, const char *a_values_vector_string,
        const std::vector<std::pair<int, VariableType>> &a_vars_vector,
        std::array<double, NUM_VARS> &a_values_array, const T a_default_value)
    {
        // how many values do I need to get?
        int num_values = a_vars_vector.size();
        // make a container for them, and load
        std::vector<T> vars_values(num_values, a_default_value);
        pp.load(a_values_vector_string, vars_values, num_values, vars_values);

        // populate the values_array for the NUM_VARS values with those read in
        a_values_array.fill(a_default_value);
        for (int i = 0; i < num_values; i++)
        {
            int icomp = a_vars_vector[i].first;
            CH_assert(a_vars_vector[i].second == VariableType::evolution);
            a_values_array[icomp] = vars_values[i];
        }
    }

    // General parameters
    int verbosity;
    std::vector<int> N_vect;
    double L;                               // Physical sidelength of the grid
    std::array<double, CH_SPACEDIM> center; // grid center
    IntVect ivN;                 // The number of grid cells in each dimension
    double coarsest_dx;          // The coarsest resolution
    int max_level;               // the max number of regriddings to do
    int num_ghosts;              // must be at least 3 for KO dissipation
    int tag_buffer_size;         // Amount the tagged region is grown by
    Vector<int> ref_ratios;      // ref ratios between levels
    Vector<int> regrid_interval; // steps between regrid at each level
    int max_steps;
    bool ignore_checkpoint_name_mismatch;   // ignore mismatch of variable names
                                            // between restart file and program
    double dt_multiplier, stop_time;        // The Courant factor and stop time
    int checkpoint_interval, plot_interval; // Steps between outputs
    int max_grid_size, block_factor;        // max and min box sizes
    double fill_ratio; // determines how fussy the regridding is about tags
    std::string checkpoint_prefix, plot_prefix; // naming of files
    bool write_plot_ghosts;
    int num_plot_vars;
    std::vector<std::pair<int, VariableType>>
        plot_vars; // vars to write to plot file

    // Boundary conditions
    std::array<bool, CH_SPACEDIM> isPeriodic;     // periodicity
    BoundaryConditions::params_t boundary_params; // set boundaries in each dir
    bool nonperiodic_boundaries_exist;
    bool boundary_solution_enforced;

    // For tagging
    double regrid_threshold;
};

#endif /* CHOMBOPARAMETERS_HPP_ */
