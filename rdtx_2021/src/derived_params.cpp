#ifndef DERCON_C_
#define DERCON_C_
#include <iostream>
#include <sstream>
#include <math.h>
#include "global_constants.h"
#include "constants.h"
#include "input_params.h"

//------------------------------------------------------------
// ****** DERIVED CONSTANTS **********
//------------------------------------------------------------

long double totalErad = 0.0;
long double phasew = 0.0;

fourtens metrictens; // The metric tensor
std::string GetstrTime(long double time);
std::string Data_directory = "output";
std::string Static_Run_directory = "/run";
std::string Run_directory = Data_directory + Static_Run_directory;
std::string Field_directory = Run_directory + "/fields/";
std::string filename_field = Field_directory + "fields_";
std::string filename_potential = Field_directory + "potentials_";
std::string filename = Run_directory + "/output_par_";
std::string filename_summary = Run_directory + "/output_all.dat";
std::string filename_spectrum = Run_directory + "/spectrum_theta_0.dat";
std::string filename_spec = Run_directory + "/spectrum_par_";

#endif
