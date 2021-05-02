/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_Ham,

    c_Mom,

    c_Weyl4_Re,

    c_Weyl4_Im,

    c_rho,

    c_proca_Re,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "Ham",
    "Mom",
    "Weyl4_Re", 
    "Weyl4_Im",
    "rho",
    "proca_Re"
    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */
