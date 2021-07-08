

#ifndef FIELDFUNC_C_
#define FIELDFUNC_C_
#include "fourvec.h"
#include "constants.h"
#include "headers.h"

namespace with_focusing
{
  //For gaussian beam, waist and Guoy phase
  long double beam_waist_sq(fourvec &x)
  {
    long double z_temp = ((*laserpulse1::dir0) * x) - laserpulse1::z_focus;
    return w0sq * (1.0L + z_temp * z_temp * invzRsq);
  }
  long double inv_beam_Radius(fourvec &x)
  {
    long double z_temp = ((*laserpulse1::dir0) * x) - laserpulse1::z_focus;
    return z_temp / (z_temp * z_temp + laserpulse1::zR * laserpulse1::zR);
  }
  long double Gouy_phase(fourvec &x)
  {
    long double z_temp = ((*laserpulse1::dir0) * x) - laserpulse1::z_focus;
    return atan2(z_temp, laserpulse1::zR);
  }

  //For gaussian beam, waist and Guoy phase
  long double beam_waist_sq2(fourvec &x)
  {
    long double z_temp = ((*laserpulse2::dir0) * x) - laserpulse2::z_focus;
    return w0sq2 * (1.0 + z_temp * z_temp * invzRsq2);
  }
  long double inv_beam_Radius2(fourvec &x)
  {
    long double z_temp = ((*laserpulse2::dir0) * x) - laserpulse2::z_focus;
    return z_temp / (z_temp * z_temp + laserpulse2::zR * laserpulse2::zR);
  }
  long double Gouy_phase2(fourvec &x)
  {
    long double z_temp = ((*laserpulse2::dir0) * x) - laserpulse2::z_focus;
    return atan2(z_temp, laserpulse2::zR);
  }
}
namespace no_focusing
{
  //For gaussian beam, waist and Guoy phase
  long double beam_waist_sq(fourvec &x)
  {
    return w0sq;
  }
  long double inv_beam_Radius(fourvec &x)
  {
    return 0.0L;
  }
  long double Gouy_phase(fourvec &x)
  {
    return 0.0L;
  }

  // For pulse 2
  //For gaussian beam, waist and Guoy phase
  long double beam_waist_sq2(fourvec &x)
  {
    return w0sq2;
  }
  long double inv_beam_Radius2(fourvec &x)
  {
    return 0.0L;
  }
  long double Gouy_phase2(fourvec &x)
  {
    return 0.0L;
  }
}

namespace asymmetric_gaussian_pulse
{
  long double a_laser_t_env(fourvec &x, long double phase_in_exponent, long double inv_t0squared)
  {

    long double exponent = 0.0L;

    exponent = -phase_in_exponent * phase_in_exponent * inv_t0squared;
    if (phase_in_exponent < 0.0L)
      exponent *= laserpulse1::pulse_front_steepness;
    long double ans = exp(exponent);
    // if (x.gett()<t0) ans=ans*0.0; //step function
    return ans;
  }
}

namespace truncated_gaussian_pulse
{
  long double a_laser_t_env(fourvec &x, long double phase_in_exponent, long double inv_t0squared)
  {

    long double exponent = 0.0L;

    exponent = -phase_in_exponent * phase_in_exponent * inv_t0squared;
    //if (phase_in_exponent<0.0) exponent*=laserpulse1::pulse_front_steepness;
    long double ans = exp(exponent);
    if (fabs(exponent) <= 3.0L)
    {
      // do nothing
    }
    else
    {
      ans = 0.0L;
    }
    // if (x.gett()<t0) ans=ans*0.0; //step function
    return ans;
  }
}

namespace square_pulse
{
  long double a_laser_t_env(fourvec &x, long double phase_in_exponent, long double inv_t0squared)
  {

    long double exponent = 0.0L;
    exponent = phase_in_exponent * phase_in_exponent * inv_t0squared;
    if (exponent <= 0.5L)
    {
      exponent = 1.0L;
    }
    else
    {
      exponent = 0.0L;
    }
    return exponent;
  }
}

// Electromagnetic Tensor Fab
// Vec potential a, position x, direction xhat
fourtens F(fourvec &x)
{
  fourvec xhat;
  fourvec diff_a, diff_au, diff_al;
  fourtens ans;
  fourtens ans1, ans2;
  ans.Null();

  for (int ii = 0; ii < 4; ++ii)
  {
    xhat.hat(ii);
    a_field(x + xhat * delta_step, diff_au);
    a_field(x - xhat * delta_step, diff_al);
    diff_a = (diff_au - diff_al);
    diff_a = diff_a * (0.5L / delta_step); // center diff
    diff_a.Minkowskii();
    xhat.Minkowskii();
    ans1.dyadic(xhat, diff_a);
    ans2.dyadic(diff_a, xhat);
    ans = ans + ans1 - ans2;
  }
  return ans;
}
#endif
