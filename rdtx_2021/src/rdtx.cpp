
#ifndef RDTX_C_
#define RDTX_C_

#include "headers.h"
//#include <wx/textctrl.h>
int save_defaults();
int load_defaults();

void operator<<(std::ofstream &str_out, fourvec vec);

int focus_particles(long double **par_x, long double **par_v);

int rdtx_run(fourvec *x0, fourvec *x1, fourvec *x2,
             fourvec *v0, fourvec *v1)
{

  long double start_time = clock();

  set_initial(); //!!!!
                 //----------------------------------------

  std::cout << "\ndtau = " << dtau0 << ", delta = " << delta_step << ", domega = " << domega << "\n";
  std::cout << "1/2pidtau = " << 0.5 / PI / dtau0 << ", 1/(2pi*dtau^3) = "
            << 1.0 / dtau0 / dtau0 / dtau0 / PI / 2.0 << ", 1/2pidelta_dtau = " << 0.5 / PI / dtau0
            << ", omega_max = " << omega_max << "\n";
  std::cout << "omega_max*dtau/2pi = " << 0.5 * omega_max * dtau0 / PI << ',';
  std::cout << "\nRadius of particle as seen by radiation / (c/omega0) = " << 0.5 * pow(sizeofbunch / densityofbunch * shape_factor_constant, 1.0L / 3.0L);
  std::cout << "\nv_group = " << laserpulse1::v_group;
  std::cout << "\nv_phase = " << laserpulse1::v_phase << "\n__________________________________\n\n";
  std::cout << "Laser pulse 1 PEAK a0 = " << a0_pulse1 << '\n';
  std::cout << "Laser pulse 2 PEAK a0 = " << a0_pulse2 << '\n';
  std::cout << "Virtual Spectrometer order = " << order_of_spectrometer << '\n';

  /*
     Check the directory, if it exists, make new one with new direcotry
    */
  int result = 0;
  std::string unix_com = "test -d " + Data_directory;
  result = system(unix_com.c_str());
  if (result)
  {
    unix_com = "mkdir " + Data_directory;
    result = system(unix_com.c_str());
  }
  Run_directory = Data_directory + Static_Run_directory;
  unix_com = "test -d " + Run_directory;
  result = system(unix_com.c_str());
  if (result)
  {
    unix_com = "mkdir " + Run_directory;
    result = system(unix_com.c_str());
    unix_com = "mkdir " + Field_directory;
    result = system(unix_com.c_str());
  }
  if (ifdumpnew)
  {
    result = 0;
    std::string unix_com = "test -d " + Run_directory;
    result = system(unix_com.c_str());
    if (!result)
    {
      int iii = 1, testdirectory = 0;
      std::ostringstream numtemp;
      std::string number;
      while (!testdirectory)
      {
        numtemp << iii;
        if (iii < 100)
          number = '0';
        if (iii < 10)
          number += '0';
        number += numtemp.str();
        unix_com = "test -d " + Run_directory + number;
        result = system(unix_com.c_str());

        if (result)
        {
          Run_directory = Run_directory + number;
          unix_com = "mkdir " + Run_directory;
          result = system(unix_com.c_str());
          Field_directory = Run_directory + "/fields/";
          unix_com = "mkdir " + Field_directory;
          result = system(unix_com.c_str());
          testdirectory = 1;
        }
        ++iii;
        numtemp.str("");
      }
    }
  }

  filename = Run_directory + "/output_par_";
  filename_summary = Run_directory + "/output_all.dat";
  filename_spectrum = Run_directory + "/spectrum_theta_0.dat";
  filename_spec = Run_directory + "/spectrum_par_";
  Field_directory = Run_directory + "/fields/";
  filename_field = Field_directory + "fields_";
  filename_potential = Field_directory + "potentials_";

  //initial four velocity equiv to energy/momentum
  fourvec vi(1.0L, 0.0L, 0.0L, 0.0L);
  fourvec xi(0.0L, 0.00L, 0.00L, taumax * 0.5L);

  // Spin four vector
  fourvec Si(0.0L, 0.0L, 0.0L, 0.0L);

  /*
    Want to print out fields to see if they are ok
  */

  // write fields
  if (if_dump_fields)
  {
    std::ofstream outfile_fields;
    std::ofstream outfile_potentials;
    long double ftime = 0.0;
    long double Nfieldstepsmone = n_field_dumps - 1;
    if (Nfieldstepsmone < 1.0)
    {
      Nfieldstepsmone = 1.0;
    }
    long double dftime = time_max / (Nfieldstepsmone);

    std::cout << "Writing fields/potentials...\n";
    for (int ii = 0; ii < n_field_dumps; ++ii)
    {
      ftime = ii * dftime;
      std::cout << "t = " << ftime << "\n";
      // int to string
      std::ostringstream numbuffer;
      numbuffer << ii;
      std::string field_string = filename_field + numbuffer.str() + ".dat";
      outfile_fields.open(field_string.c_str());
      std::string potential_string = filename_potential + numbuffer.str() + ".dat";
      outfile_potentials.open(potential_string.c_str());

      // now dump fields and potentials
      dump_fields(outfile_fields, outfile_potentials, ftime);
      outfile_fields.close();
      outfile_potentials.close();
    }
  }

  calculatedamping(); // gaussian filter for spectrum

  /*
    Initiate particle positions and momenta
   */
  if (ionization_on)
  { // then starts at rest!
    bunch_momentum[0] = 0.0L;
    bunch_momentum[1] = 0.0L;
    bunch_momentum[2] = 0.0L;
  }
  initialize_particles(par_x, par_v, par_S);
  // add particle focusing
  if (parzfocus)
  {
    focus_particles(par_x, par_v);
  }
  /*
    Now run for all particles
   */
  std::ofstream outfile_summary;
  outfile_summary.open(filename_summary.c_str());
  outfile_summary << "Summary of all particle data\n0\n"
                  << Npar << "\n";
  outfile_summary << "par#\tt0\tx0\ty0\tz0\tE0\tpx0\tpy0\tpz0\t";
  outfile_summary << "tf\txf\tyf\tzf\tEf\tpxf\tpyf\tpzf\tERad\tPmax/Plas\n";

  long double final_time = 0.0L;
  ResetBins();
  ResetSpec();
  spectrumerror = 0; // if there is an error in spectral calculation - returns from SpectralCalculation-only once
  final_time = 0.0L;
  totalErad = 0.0;
  for (int i = 0; i < Npar; ++i)
  {

    if (!if_import_data)
    {
      std::cout << "PARTICLE " << i + 1 << " of " << Npar << "\n";

      xi.sett(par_x[0][i]);
      vi.sett(par_v[0][i]);
      Si.sett(par_S[0][i]);
      xi.setx(par_x[1][i]);
      vi.setx(par_v[1][i]);
      Si.setx(par_S[1][i]);
      xi.sety(par_x[2][i]);
      vi.sety(par_v[2][i]);
      Si.sety(par_S[2][i]);
      xi.setz(par_x[3][i]);
      vi.setz(par_v[3][i]);
      Si.setz(par_S[3][i]);

      std::cout << "STEP 1: Integrating equations of motion for particles\n\n";

      // Reset ionizatoin state
      particleisfree = 1.0L * !ionization_on;

      // Calculate electron trajectories

      final_time = par_traj(vi, xi, Si, i, outfile_summary);
    }

    if (if_calc_spectrum || (if_import_data))
    {
      std::cout << "STEP 2: Calculating angular distribution of radiation\n";

      theta_scatter = theta_scatter_min; // very important!!!
      for (int thetabin = 0; thetabin < N_theta_bins; ++thetabin)
      {

        std::cout << "Theta = " << theta_scatter << ", " << thetabin << " of " << N_theta_bins << "\n";

        sinthetascatter = sin(theta_scatter);
        costhetascatter = cos(theta_scatter);
        excur_y = costhetascatter * costhetascatter;
        excur_z = sinthetascatter * sinthetascatter;
        fourvec Kspec0(1.0L, 0.0L, sinthetascatter, costhetascatter);
        *Kspec = Kspec0;

        if (!if_import_data)
        {
          SpectralEndpoint(thetabin, 0.0, int_v[0], int_x[0], 1.0);
          if (SpectralCalculation(thetabin, i))
          {
            return 0;
          } // go back to gui
          SpectralEndpoint(thetabin, taumax, int_v_m[Nmax], int_x[0], -1.0);
        }
        else
        {
          SpectralEndpoint(thetabin, 0.0, v0[0], x0[0], 1.0);
          SpectralEndpoint(thetabin, taumax, v0[Nmax], x0[Nmax], -1.0);
        }
        theta_scatter += dtheta;
      }
      theta_scatter = theta_scatter_min;
      for (int thetabin = 0; thetabin < N_theta_bins; ++thetabin)
      {
        //      if (usercancelrun(control,outstrstream,"Calculating spectrum")){return 0;}
        if (!coherent_addition)
        {
          Calc_Spec(thetabin);
        }
        theta_scatter += dtheta;
      }
      if (!coherent_addition)
      {
        ResetBins();
      }
    }
  }

  outfile_summary.close();

  if (if_calc_spectrum || (if_import_data))
  {
    std::cout << "Final combination of spectral components...\n";

    theta_scatter = theta_scatter_min;
    for (int thetabin = 0; thetabin < N_theta_bins; ++thetabin)
    {
      if (coherent_addition)
      {
        Calc_Spec(thetabin);
      }
      theta_scatter += dtheta;
    }

    // Normalize spectrum correctly.
    theta_scatter = theta_scatter_min;
    for (int thetabin = 0; thetabin < N_theta_bins; ++thetabin)
    {

      for (int i = 0; i < N_omega_bins; ++i)
      {
        SpecInt_para[i][thetabin] *= particle_shape_k[i] * omega[i] * omega[i] * spectrum_const;
        SpecInt_perp[i][thetabin] *= particle_shape_k[i] * omega[i] * omega[i] * spectrum_const;
      }
      theta_scatter += dtheta;
    }
    theta_scatter = theta_scatter_min;
    for (int thetabin = 0; thetabin < N_theta_bins; ++thetabin)
    {

      std::cout << "Theta = " << theta_scatter << ", " << thetabin << " of " << N_theta_bins << "\n";

      sinthetascatter = sin(theta_scatter);
      costhetascatter = cos(theta_scatter);
      excur_y = costhetascatter * costhetascatter;
      excur_z = sinthetascatter * sinthetascatter;

      Write_Spec(-1, thetabin);
      theta_scatter += dtheta;
    }
  }
  std::cout << "finished\n";
  long double end_time = clock();
  std::cout << "Total radiated energy = " << totalErad * 511.0 << " keV\n";
  std::cout << "\nEND OF CALCULATION\n";
  std::cout << "Calculation took " << GetstrTime(end_time - start_time) << '\n';

  return 0;
}

#endif
