/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "ChomboParameters.hpp"
#include "SimulationParametersBase.hpp"
#include "GRParmParse.hpp"
// Problem specific includes:
#include "KerrSchildFixedBG.hpp"
#include "Potential.hpp"
#include "SpheroidalExtraction.hpp"
#include "InitialConditions.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

inline void read_number(double &data, std::string path){

    double x;
    std::ifstream inFile;
    std::string filename = "check";

    inFile.open(path);
    if (!inFile) {
        std::cout << "Unable to open file";
        exit(1); // terminate with error
    }
    inFile >> data;

    inFile.close();


}

inline void read_file(std::vector<double> &data, std::string path, bool verbose = false  ){

    double x;
    std::ifstream inFile;
    std::string filename = "check";

    inFile.open(path);
    if (!inFile) {
        std::cout << "Unable to open file";
        exit(1); // terminate with error
    }

    while (inFile >> x) {
        data.push_back(x);
    }

    inFile.close();

    if(verbose){
        for (std::vector<double>::iterator it = data.begin(); it != data.end(); ++it) {
            std::cout << *it << " "<<std::endl;
        }
    }
}


class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("nan_check", nan_check, 1);
        pp.load("sigma", sigma, 0.1);
        pp.load("integral_filename", integral_filename);

        // Initial and Kerr data
        pp.load("bh_mass", bg_params.mass);
        pp.load("bh_spin", bg_params.spin);
        pp.load("bh_center", bg_params.center, center);
        pp.load("proca_mass", potential_params.mass);
        pp.load("proca_self_interaction", potential_params.self_interaction);
        pp.load("proca_damping", proca_damping);
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("threshold_phi",threshold_phi);
        pp.load("threshold_K",threshold_K);


        field_mu = potential_params.mass;

        // extraction params
        dx.fill(coarsest_dx);
        origin.fill(coarsest_dx / 2.0);

        std::string folder = "vector_field_star_self_interacting/Lambda_0/Dim_4/f0_0.165/cA4_0.0/";
        std::string filename_a0 = "a0.dat";
        std::string filename_da0dr = "da0dr.dat";
        std::string filename_a1 = "a1.dat";
        std::string filename_m = "m.dat";
        std::string filename_sig = "sigma.dat";
        std::string filename_rvals = "rvals.dat";
        std::string filename_omega = "omega.dat";

        std::string path_a0 = folder + filename_a0;
        std::string path_da0dr = folder + filename_da0dr;
        std::string path_a1 = folder + filename_a1;
        std::string path_m  = folder + filename_m;
        std::string path_sig  = folder + filename_sig;
        std::string path_rvals  = folder + filename_rvals;
        std::string path_omega  = folder + filename_omega;

        read_file(initalcondition_data.a0,path_a0);
        read_file(initalcondition_data.da0dr,path_da0dr);
        read_file(initalcondition_data.a1,path_a1);
        read_file(initalcondition_data.m,path_m);
        read_file(initalcondition_data.sig,path_sig);
        read_file(rvals,path_rvals);
        read_number(initalcondition_data.omega,path_omega);

        pp.load("centerProca1", initalcondition_data.center1);
        pp.load("centerProca2", initalcondition_data.center2);

        initalcondition_data.spacing = rvals[1] - rvals[0];

        // Extraction params
        pp.load("num_extraction_radii", extraction_params.num_extraction_radii,
                1);
        // Check for multiple extraction radii, otherwise load single
        // radius/level (for backwards compatibility).
        if (pp.contains("extraction_levels"))
        {
            pp.load("extraction_levels", extraction_params.extraction_levels,
                    extraction_params.num_extraction_radii);
        }
        else
        {
            pp.load("extraction_level", extraction_params.extraction_levels, 1,
                    0);
        }
        if (pp.contains("extraction_radii"))
        {
            pp.load("extraction_radii", extraction_params.extraction_radii,
                    extraction_params.num_extraction_radii);
        }
        else
        {
            pp.load("extraction_radius", extraction_params.extraction_radii, 1,
                    0.1);
        }
        pp.load("num_points_phi", extraction_params.num_points_phi, 2);
        pp.load("num_points_t", extraction_params.num_points_t, 5);
        if (extraction_params.num_points_t % 2 == 0)
        {
            extraction_params.num_points_t += 1;
            pout() << "Parameter: num_points_t incompatible with Simpson's "
                   << "rule so increased by 1.\n";
        }
        pp.load("extraction_center", extraction_params.center, center);
        pp.load("zaxis_over_xaxis", extraction_params.zaxis_over_xaxis, 1.0);
        pp.load("write_extraction", extraction_params.write_extraction, false);
    }

    // Problem specific parameters
    double field_mu;
    double sigma, proca_damping;
    double threshold_phi, threshold_K;
    int nan_check;
    std::array<double, CH_SPACEDIM> origin,
        dx; // location of coarsest origin and dx
    std::vector<double> rvals;
    InitialConditions::params_t initalcondition_data;
    double G_Newton;
    double spacing;
    double omega;
    std::string integral_filename;
    // Collection of parameters necessary for the sims
    KerrSchildFixedBG::params_t bg_params;
    Potential::params_t potential_params;
    SpheroidalExtraction::params_t extraction_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
