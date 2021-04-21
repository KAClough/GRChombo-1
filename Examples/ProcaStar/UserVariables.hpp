/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "DiagnosticVariables.hpp"
#include "CCZ4UserVariables.hpp"

// assign an enum to each variable
enum
{
    c_Avec0_Re = NUM_CCZ4_VARS,

    c_Avec1_Re,
    c_Avec2_Re,
    c_Avec3_Re,

    c_Evec1_Re,
    c_Evec2_Re,
    c_Evec3_Re,

    c_Zvec_Re,

    c_Avec0_Im,

    c_Avec1_Im,
    c_Avec2_Im,
    c_Avec3_Im,

    c_Evec1_Im,
    c_Evec2_Im,
    c_Evec3_Im,

    c_Zvec_Im,

    NUM_VARS
};

namespace UserVariables
{
    static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS> user_variable_names = {
    "phi_Re",

    "Avec1_Re", "Avec2_Re", "Avec3_Re",

    "Evec1_Re", "Evec2_Re", "Evec3_Re",

    "Z_Re",

    "phi_Im",

    "Avec1_Im", "Avec2_Im", "Avec3_Im",

    "Evec1_Im", "Evec2_Im", "Evec3_Im",

    "Z_Im"

    };

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
