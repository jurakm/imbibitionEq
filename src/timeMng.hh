/**
 * timeMng.hh
 *
 *  Class for management of time increments.
 */

#ifndef SRC_TIMEMNG_HH_
#define SRC_TIMEMNG_HH_

#include "parameters.hh"

/// \brief  Time manager class.
///
/// Advance time by given time step respecting the end time and the output times.
/// Change the time step depending on the number of solver iterations.

template <typename Real>
struct TimeMng{
  /** Current time. */
  Real time;
  /** Current time increment. */
  Real dt;
  /** Requested time increment for the next step. */
  Real requested_dt;
  /** Maximum time increment. */
  Real dtmax;
  /** Time increment used for writing the output files. */
  Real dtout;
  /** Time of the next output operation. */
  Real output_time;
  /** End time of simulation. */
  Real tend;
  /** Count of time steps. */
  int  count;
  /** Count of output operations (start with zero). */
  int  output_count;
  /** Tolerance for real number comparison. */
  const Real TOL = 1.0e-10;
  /**
   * Initialize the time manager from the parameters object.
   */
  template <typename Parameters>
  TimeMng(Parameters const & params){
    time = 0.0;
    output_time = 0.0;
    dt    = params.dt;
    requested_dt = dt;
    dtmax = std::max(dt, params.dtmax);
    dtout = std::max(dt,params.dtout);
    tend  = std::max(dtout, params.tend);
    count = 0;
    output_count = 0;
  }

  /**
   * Manage requested time step as a function of a number of nonlinear iterations.
   */
  void set_requested_dt(int iterCount){
    if(iterCount < 4) requested_dt *= 1.2;
    if(iterCount > 8) requested_dt /= 1.4;
    if(requested_dt > dtmax) requested_dt = dtmax;
  }

  /**
   * Increment time for requested time step. Do not step over output or final times.
   */
  TimeMng & operator++(){
    if(time >= tend) return *this;

    ++count;
    // At this instance the time is old time. If output time
    // is reached, increment it!
    if(time >= output_time) {
	// increment output time but do not step over final time
	Real new_output_time = output_time + dtout;
	if(new_output_time > tend){
	    dtout = tend - output_time;
	    output_time = tend;
	}
	else
	    output_time = new_output_time;  // regular case

        ++output_count;
    }
    // Increment time but do not step over output time.
    // In that way the final time will also be respected.

    dt = requested_dt; // when incrementing actually used dt in the last step is not important

    Real new_time = time + dt;
    if(new_time > output_time) {
	dt = output_time - time;
	time = output_time;
    }
    else if(new_time + 0.1*dt >= output_time){
	// Here we increment dt (less than 10 %) in order to avoid too small steps
	dt = output_time - time;
	time = output_time;
    }
    else
      time = new_time;   // regular case

    return *this;
  }

  /** Postfix form of the increment. Same as the prefix form. */
  TimeMng & operator++(int){
    ++*this;
    return *this;
  }

  /** Is it time for output? */
  bool doOutput() const {
    if ( std::abs(time - output_time) < TOL ) return true;
    return false;
  }

  /** Is the final time reached? */
  bool done() const {
    if ( std::abs(time - tend) < TOL || time > tend) return true;
    return false;
  }


};



#endif /* SRC_TIMEMNG_HH_ */
