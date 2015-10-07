/*
 * linear_analytic.hh
 *
 *  Created on: 25 Jun 2015
 *      Author: jurak
 */

#ifndef SRC_LINEAR_ANALYTIC_HH_
#define SRC_LINEAR_ANALYTIC_HH_

#include <iostream>
#include <cmath>
#include <string>
#include <functional>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>

#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/exception.h>

#include "parameters.hh"
#include "mgf.hh"
#include "tableXY.hh"

// Find interval (0,M) on which the 1-erf() function falls under
// given  tolerance. Practically the result is M=6.
double find_upper_bound(double TOL){
  double M = 1.0;
  for(int i=0; i < 30; ++i){
    if(std::abs(1.0 - std::erf(M)) > TOL) M= M+1.0;
    else break;
  }
  return M;
}

// Calculate y=Y(x) by linear interpolation.
//double interpolate(boost::numeric::ublas::vector<double> const & X,
//		           boost::numeric::ublas::vector<double> const & Y, double x,
//		           unsigned int last_index) {
//	const double TOL = 1E-8;
//	double x_max = X[last_index];
//	if (x > x_max + TOL) {
//		std::cerr << "x= " << x << ", x_max = " << x_max << ", idx = " << last_index << std::endl;
//		throw std::runtime_error("interpolate: value out of bounds, error 541.");
//	}
//	// Small negative numbers are ossible.
////	if(x < 0.0 and -x < TOL) x = 0.0;
//	auto it = std::upper_bound(X.begin(), X.begin() + last_index + 1, x);
//	auto index = it - X.begin(); // this is upper bound
//	if (index == 0){
////		std::cout << "x = " << x << ", x_max = " << x_max << ", idx = " << last_index << std::endl;
////		int bbb; std::cin >> bbb;
////		return x;
//		throw std::logic_error("Internal error 542.");
//	}
//	if (index > last_index) {
//		// There is no bigger element, return largest value.
//		// This can happen only on the last element.
//		assert(index == last_index + 1);
//		return Y[index - 1];
//	}
//	const double y1 = Y[index - 1];
//	const double y2 = Y[index];
//	const double x1 = X[index - 1];
//	const double x2 = X[index];
//	double y = y1;
//	if (std::abs(x2 - x1) >= TOL)
//		y += (x - x1) * (y2 - y1) / (x2 - x1);
//	return y;
//}






/** \brief Class for calculation of analytic solution in 1D.
 *
 * Class holding the integrand for calculation of the integral
 * (boundary layer function)
 *   \f[
 *   Z(\xi,t) = (2/\sqrt{\pi}) \int_{x/(2\delta\sqrt{\tau(t)})}^\infty ((g\circ t)(\tau(t)-\xi^2/(4 v^2)) -g(0))exp(-v^2) dv.
 *   \f]
 *   where \f$\xi\f$ is  a boundary value variable, that is, at \f$x=0\f$
 *   \f[
 *   \xi = x/\overline{\delta},\quad \overline{\delta} = \delta \sqrt{\frac{k_m\overline{\alpha}_m}{\Phi_m}}.
 *   \f]
 * For non wetting phase flux \f$Q_n\f$ we have the formula
 * \f[
 *     Q_n = \frac{4}{\sqrt{\pi}}\overline{\delta}\Phi_m a(t)
 *                         \int_0^{\sqrt{\tau(t)}} (g\circ t)'(\tau(t)-\lambda^2)\, d\lambda.
 * \f]
 * Different cases:
 *   Constant approximation: \f$ a(t) = 1,\; \tau(t) = t. \f$
 *   Variable approximation: \f$ a(t) = \alpha(g(t))/\alpha_m,\; \tau(t) = \int_0^t a(u)du. \f$
 */
template <typename Params>
class AnalyticSolution{
  public:
	/// Constructor.
    AnalyticSolution(Params const & params);
    /** Calculate the flux in the case of constant linearisation. The result is
     *  given on an equidistant time mesh and is stored in lin_flux variable.
     */
    void calculate_flux();
    // beta version
    void calculate_flux_beta();
    /** Print the flux calculated in calculate_linear_const_flux() on a given stream.
     *  If the flux is not already calculated, the function will call calculate_linear_const_flux().
     */
    void print_flux(std::ostream & out);
    void print_flux_beta(std::ostream & out);
    /** Calculate the solution of constant linearization in the time instant time.
     *  The solution is stored in lin_const_solution variable.
     */
    void calculate_solution(double time);
    /** Print the solution of constant linearization calculated by calculate_linear_const_solution()
     *  method. We do not verify that the solution is indeed calculated.
     */
    void print_solution(std::ostream & out);
     /** Print tau(t) in file for debugging. */
    void print_tau(std::ostream & out);
//    void print_tau1(std::ostream & out);
    /** tau( t ). */
    double tau(double t) const;
    /** t= t(tau), invers of tau. */
    double inv_tau(double tau) const;
    /** Model we are simulating. */
    int model() const { return model_; }
  private:
    const Params & params_;
    const int model_;   // capture params_.model since it can change during the simulation

    void set_x(double xx) { x_=xx;  xi_ = x_/scaled_delta; }

    /** Integrand in calculation of BL function Z(x,t) in linear constant approximation case. */
    double g_shifted(double v, double t) const {
      return factor * ( (bdry(tau(t) - xi_*xi_/(4*v*v)) - bdry(0.0)) * std::exp(-v*v) );
    }
    /** Integrand in calculation of the nonwetting phase flux, scalled with the right factor.
     * Linear constant case. */
    double dg_dt_shifted(double v, double t){ return factor_der *  dg_dt(tau(t) - v*v); }
    double dbetag_dt_shifted(double v, double t){ return factor_der_beta *  dbetag_dt(tau(t) - v*v); }
    /** Integral lower bound. */
    double lower_bound(double t) const { return xi_/(2*std::sqrt(tau(t)));}
    /** Boundary condition at the time t. */
    double bdry(double t) const {
    	double val = 0.0;
    	if(model_ == Params::analytic_const)	val = params_.bdry(t);
    	else                                    val = params_.bdry(inv_tau(t)); // t is tau here
//    	std::cout << t << " " << val << "\n";
    	return val;
    }
    // This is necessary for the second form of the flux
    double beta_bdry(double t) const {
        double val = 0.0;
        if(model_ == Params::analytic_const)	val = params_.beta(params_.bdry(t));
        else                                    val = params_.beta(params_.bdry(inv_tau(t))); // t is tau here
    //   	std::cout << t << " " << val << "\n";
        return val;
    }
    /** Scaled composition alpha( g(t) ). */
    double a_g(double t) const {
    	static const double TOL = 1.0e-5;
    	//static const double bs0 = params_.beta( params_.bdry(0.0) );
    	//  najbolji rezultati su sa 0.75.
    	double val = 1.0;
    	if(model_ == Params::analytic_const)
    		val = 1.0;
    	else if(model_ == Params::analytic_var)
    	   val = params_.alpha( params_.bdry(t) )/params_.mean_alpha;
    	else if(model_ == Params::analytic_new){
    		// 0.5 daje preveliki flux, 1 daje premali.
    		const double Yt = params_.bdry(t);
    		const double dS = Yt - params_.bdry(0.0);
    		if(std::abs(dS*theta) > TOL){
    	        val = ( params_.beta( Yt ) - params_.beta( Yt - dS*theta) )/ (dS*theta);
    		}
    		else
    			val = params_.alpha(Yt);

    		val /= params_.mean_alpha;
    	}
    	else if(model_ == Params::analytic_new1){
    		const double Yt = params_.bdry(t);
    		const double Yt0 = (t>dt_bdry) ? params_.bdry(t - dt_bdry) : params_.bdry(0.0);
    		const double dS = Yt - Yt0;
    		if(std::abs(dS) > TOL){
    	        val = ( params_.beta( Yt ) - params_.beta( Yt0) )/ dS;
    		}
    		else
    			val = params_.alpha(Yt);

    		val /= params_.mean_alpha;
    	}
    	return val;
    }
    /** Derivative of the boundary value function.  */
    double dg_dt(double t){
    	double val = 0.0;
    	if(t < h) val = (bdry(t+h) - bdry(t))/h;
    	else val = (bdry(t+h) - bdry(t-h))/(2*h);
    	return val;
    }

    /** Derivative of the boundary value function -- beta version.  */
       double dbetag_dt(double t){
       	double val = 0.0;
       	if(t < h) val = (beta_bdry(t+h) - beta_bdry(t))/h;
       	else val = (beta_bdry(t+h) - beta_bdry(t-h))/(2*h);
       	return val;
       }

    /** Calculate new entry in  new_time_field by calculating
     * \f[ \int_{t_0}^t \alpha_m( g(s)) ds,\f]
     * where \f$t_0\f$ is given by last_time_index_.
     */
    void integrate_alpha_bdry(double t);

    double scaled_delta = 0.0;
    double x_  = 0.0;
    double xi_ = 0.0;
    const double factor = 2.0/std::sqrt(M_PI);
    double factor_der = 0.0;
    double factor_der_beta = 0.0;
    const double h = 1.0E-7;
    double dt_table_ = 0.0;
    double theta = 1.0;
    double dt_bdry = 0.0;

    /** Linear model flux; pairs (t, flux(t)).  */
    std::vector<std::pair<double,double>> lin_flux;
    // beta version
    std::vector<std::pair<double,double>> lin_flux_beta;
    /** Linear constant solution; pairs (x, S_w(x,t)) at fixed time t.  */
    std::vector<std::pair<double,double>> lin_solution;
    /** Last index of calculated time in new_time_field. */
//    unsigned int last_time_index_ = 0;
//    const unsigned int default_table_size_ = 10000;
    /** time_ holds time  t. */
//    boost::numeric::ublas::vector<double> time_;
    /** tau_time_ holds tau(t): tau(t) = int_0^t alpha(g(s))ds.   */
//    boost::numeric::ublas::vector<double> tau_time_;

    /** x-coordinates for the solution. */
	std::vector<double> x_coo;

	TableXY<boost::numeric::ublas::vector<double> > tau_table;
};

template<typename Params>
AnalyticSolution<Params>::AnalyticSolution(Params const & params) :
		params_(params), model_(params_.model) {
	// take parameters
	double delta = params.delta;
	double perm = params.k;
	double poro = params.poro;
	double mean_alpha = params.mean_alpha;
	theta = params.theta;
	dt_bdry = params.dt_bdry;
	scaled_delta = params.scaled_delta;
//   	scaled_delta = delta * std::sqrt(perm*mean_alpha/poro);
	assert(scaled_delta > 0.0);
	// factor = 2.0/std::sqrt(M_PI);
	factor_der = 2 * poro * scaled_delta * factor;
	factor_der_beta = 2 * factor * delta * std::sqrt(perm*poro/mean_alpha);

	double tend = params.tend;
	// live 10 % of space since a number of points generated is not exactly default_table_size
	dt_table_ = tend / (0.9 * tau_table.size());
//	last_time_index_ = 0;
//
//	time_.resize(default_table_size_);
//	tau_time_.resize(default_table_size_);
//
//	time_(0) = 0.0;
//	tau_time_(0) = 0.0;

	tau_table.push_back(0.0,0.0);

	// fill the tables -- calculates tau(t)
	integrate_alpha_bdry(tend);
	MGF<Params> d1(params_);
	d1.double_side_interval(x_coo);
}



template <typename Params>
double AnalyticSolution<Params>::tau(double t) const{
	// In analytic_const mode tau(t) = t and invers function t=t(tau) is simply identity.
	double val = 0.0;
//	if(model_ == Params::analytic_const)	val = t;
//	else val = interpolate(time_, tau_time_, t, last_time_index_);

	if(model_ == Params::analytic_const)	val = t;
	else val = tau_table.interpolate_y(t);

	return val;
}

template <typename Params>
double AnalyticSolution<Params>::inv_tau(double tau) const{
	// In analytic_const mode tau(t) = t and invers function t=t(tau) is simply identity.
	double val = 0.0;
//	if(model_ == Params::analytic_const)	val =  tau;
//	else val = interpolate(tau_time_, time_,  tau, last_time_index_);

	if(model_ == Params::analytic_const)	val = tau;
	else val = tau_table.interpolate_x(tau);

	return val;
}

template <typename Params>
void AnalyticSolution<Params>::integrate_alpha_bdry(double tend){
//	unsigned int old_last_index_time = last_time_index_;
//	double t_inf =  time_[last_time_index_];
//	double integr=  tau_time_[last_time_index_];


	double t_inf  = tau_table.get_last_x();
	double integr = tau_table.get_last_y();

	if(t_inf >= tend){
		std::cerr << "Integrand<Params>::integrate_alpha_bdry(double t)  -- t_inf >= t!\n";
		// there is no need for integration -- this situation is possibly an error.
		return;
	}
//	auto tmp = o2scl::err_hnd; -- do not return old handler
	auto new_err_hnd = new o2scl::err_hnd_cpp();
	o2scl::err_hnd =  new_err_hnd;
	o2scl::funct11 f = std::bind(std::mem_fn<double(double)const>(&AnalyticSolution<Params>::a_g),
	    		                 this, std::placeholders::_1);

    o2scl::inte_adapt_cern<o2scl::funct11, 2000> inte_formula;
//	o2scl::inte_qag_gsl<> inte_formula;
	inte_formula.tol_abs = inte_formula.tol_rel = 1e-5;
	double res=0.0, err=0.0;
	double intpart;
	double fractpart = std::modf((tend - t_inf)/dt_table_, &intpart);
	if(fractpart >= 0.3) intpart++;

	for(unsigned int i=0; i<intpart; ++i){
		double time = t_inf + (i+1)*(tend-t_inf)/intpart;
      	inte_formula.integ_err(f,t_inf,time,res,err);
//     	last_time_index_++;
//     	if(last_time_index_ >= time_.size()){
//     		time_.resize(2*last_time_index_);
//     		tau_time_.resize(2*last_time_index_);
//     	}
//      	time_[last_time_index_] = time;
//      	tau_time_[last_time_index_] = integr + res;

      	tau_table.push_back(time, integr + res);
	}
    if( !tau_table.check_order() )
    	throw std::runtime_error("AnalyticSolution<Params>::integrate_alpha_bdry: ordering incorrect!");
//	o2scl::err_hnd = tmp;  -- no need. Poses a problem in parallel
//	delete new_err_hnd;
}

// just for debugging
//template <typename Params>
//void AnalyticSolution<Params>::print_tau(std::ostream & out){
//	out << "#  time    time    tau(time)      tau(time)\n";
//	for(unsigned int i = 0; i <= last_time_index_; ++i){
//		out << time_[i] << " " << inv_tau(tau_time_[i]) << " " << tau_time_[i] << "  "<< tau(time_[i]) << "\n";
//	}
//}
// just for debugging
template <typename Params>
void AnalyticSolution<Params>::print_tau(std::ostream & out){
	out << "#  time    time    tau(time)      tau(time)\n";
	for(int i = 0; i <= tau_table.last(); ++i){
		out << tau_table.get_x(i)  << " " << inv_tau(tau_table.get_y(i)) << " " << tau_table.get_y(i)
				<< "  "<< tau(tau_table.get_x(i)) << "\n";
	}
}

template <typename Params>
void AnalyticSolution<Params>::calculate_flux(){
	  double dt   = params_.dtout;
      int    Nsteps = params_.tend / dt;
      double t = dt;
      lin_flux.resize(Nsteps+1);
      std::fill(lin_flux.begin(), lin_flux.end(), std::make_pair(0.0,0.0));

      o2scl::funct11 fprim = std::bind(std::mem_fn<double(double,double)>(&AnalyticSolution<Params>::dg_dt_shifted),
    		                           this, std::placeholders::_1, std::cref(t));
//     auto tmp = o2scl::err_hnd; -- we do not return old handler
      auto new_err_hnd =  new o2scl::err_hnd_cpp();
      o2scl::err_hnd = new_err_hnd;
//      o2scl::inte_qag_gsl<> inte_formula;
      o2scl::inte_adapt_cern<o2scl::funct11, 2000> inte_formula;
      inte_formula.tol_rel = 1.0e-5;
      inte_formula.tol_abs = 1.0e-5;
      // Output flux -- one value for each time instant.
      // Time loop
      lin_flux[0].first = 0.0;
      lin_flux[0].second = 0.0;
      double res=0.0, err=0.0;
      for(int it=1; it <=Nsteps; ++it){
    	  res =-1; err = -1;
         // Calculate the flux by given formula
//    	 std::cout << "it = " << it << ", t = " << t << ", tau(t) = " << tau(t) << ", sqrt(tau(t))= " << std::sqrt(tau(t)) << std::endl;;
         inte_formula.integ_err(fprim,0.0,std::sqrt(tau(t)),res,err);
//    	 std::cout << ", res = " <<  res << ", err = " << err  << std::endl;
         lin_flux[it].first = t;
         lin_flux[it].second = res * a_g(t);
    //    std::cout << "err = " << err << "\n";
         t += dt;
      }

//      o2scl::err_hnd = tmp;
//      delete new_err_hnd;

      return;
}

template <typename Params>
void AnalyticSolution<Params>::calculate_flux_beta(){
	  double dt   = params_.dtout;
      int    Nsteps = params_.tend / dt;
      double t = dt;
      lin_flux_beta.resize(Nsteps+1);
      std::fill(lin_flux_beta.begin(), lin_flux_beta.end(), std::make_pair(0.0,0.0));

      o2scl::funct11 fprim = std::bind(std::mem_fn<double(double,double)>(&AnalyticSolution<Params>::dbetag_dt_shifted),
    		                           this, std::placeholders::_1, std::cref(t));
//     auto tmp = o2scl::err_hnd; -- we do not return old handler
      auto new_err_hnd =  new o2scl::err_hnd_cpp();
      o2scl::err_hnd = new_err_hnd;
//      o2scl::inte_qag_gsl<> inte_formula;
      o2scl::inte_adapt_cern<o2scl::funct11, 2000> inte_formula;
      inte_formula.tol_rel = 1.0e-5;
      inte_formula.tol_abs = 1.0e-5;
      // Output flux -- one value for each time instant.
      // Time loop
      lin_flux_beta[0].first = 0.0;
      lin_flux_beta[0].second = 0.0;
      double res=0.0, err=0.0;
      for(int it=1; it <=Nsteps; ++it){
    	  res =-1; err = -1;
         // Calculate the flux by given formula
//    	 std::cout << "it = " << it << ", t = " << t << ", tau(t) = " << tau(t) << ", sqrt(tau(t))= " << std::sqrt(tau(t)) << std::endl;;
         inte_formula.integ_err(fprim,0.0,std::sqrt(tau(t)),res,err);
//    	 std::cout << ", res = " <<  res << ", err = " << err  << std::endl;
         lin_flux_beta[it].first = t;
         lin_flux_beta[it].second = res;
    //    std::cout << "err = " << err << "\n";
         t += dt;
      }

//      o2scl::err_hnd = tmp;
//      delete new_err_hnd;

      return;
}




template <typename Params>
void AnalyticSolution<Params>::print_flux(std::ostream & out){
	// Note. Flux must be calculated first.
    for(unsigned int i = 0; i < lin_flux.size(); ++i)
    	out << lin_flux[i].first << "  " << lin_flux[i].second <<"\n";
}

template <typename Params>
void AnalyticSolution<Params>::print_flux_beta(std::ostream & out){
	// Note. Flux must be calculated first.
    for(unsigned int i = 0; i < lin_flux_beta.size(); ++i)
    	out << lin_flux_beta[i].first << "  " << lin_flux_beta[i].second <<"\n";
}

template<typename Params>
void AnalyticSolution<Params>::calculate_solution(double time) {
	double L = params_.L;
	lin_solution.resize(params_.N);

	o2scl::funct11 f = std::bind(
			std::mem_fn<double(double, double)const>(&AnalyticSolution<Params>::g_shifted),
			this, std::placeholders::_1, std::cref(time));
//	  o2scl::inte_qag_gsl<> inte_formula;
    o2scl::inte_adapt_cern<o2scl::funct11, 2000> inte_formula;
//    auto tmp = o2scl::err_hnd;
    auto new_err_hnd = new o2scl::err_hnd_cpp();
    o2scl::err_hnd =  new_err_hnd;

	double lb, ub, res1, res2, err;
	if (time == 0.0)
		for (int i = 0; i < params_.N; ++i) {
			lin_solution[i].first = x_coo[i];
			lin_solution[i].second = bdry(0.0);
		}
	else
		for (int i = 0; i < params_.N; ++i) {
			set_x(x_coo[i]);
			lb = lower_bound(time);  // depends on x!
			assert(lb >= 0.0);
			ub = std::max(6.0, lb + 1);
			// calculate Z(x,t)
			inte_formula.integ_err(f, lb, ub, res1, err);
//	        std::cout << "err = " << err <<"\n";
			set_x(L - x_coo[i]);
			lb = lower_bound(time);  // depends on x!
			ub = std::max(6.0, lb + 1);
			inte_formula.integ_err(f, lb, ub, res2, err);

			lin_solution[i].first = x_coo[i];
			lin_solution[i].second = bdry(0.0) + res1 + res2;
		}

//    o2scl::err_hnd = tmp;
//    delete new_err_hnd;
	return;
}


template <typename Params>
void AnalyticSolution<Params>::print_solution(std::ostream & out){
	for(unsigned int i = 0; i < lin_solution.size(); ++i)
	    	out << lin_solution[i].first << "  " << lin_solution[i].second <<"\n";
}

// compilation
//  g++ -std=c++11 -g linear_analytic.cc  -o analytic /usr/local/lib/libo2scl.a -lgsl -lblas
template<typename Params>
int lin_analytic_driver(Params params) {
	auto start = std::chrono::system_clock::now();
	AnalyticSolution<Params> a(params);

	double dt = params.dtout;
	int Nsteps = params.tend / dt;
	double t = 0.0;
    // 1. Calculate flux
	const std::string base( params.str_sname + params.simulation_names[a.model()] );
    const std::string tau = base + "-tau.txt";
    std::ofstream out_tau(tau);
    a.print_tau(out_tau);
    out_tau.close();
//    {
//    	    const std::string base( params.str_sname + params.simulation_names[a.model()] );
//    	    const std::string tau = base + "-tau1.txt";
//    	    std::ofstream out_tau(tau);
//    	    a.print_tau1(out_tau);
//    	    out_tau.close();
//    }

	a.calculate_flux();
	// Output flux -- one value for each time instant.
	const std::string flux( base + "-flux.txt");
	std::ofstream out_flux(flux);
	a.print_flux(out_flux);
	out_flux.close();
	a.calculate_flux_beta();
	const std::string flux_beta( base + "-flux-beta.txt");
    std::ofstream out_flux_beta(flux);
	a.print_flux_beta(out_flux_beta);
	out_flux_beta.close();

	// 2. Calculate solution
	// Time loop
	for (int it = 0; it <= Nsteps; ++it) {
		a.calculate_solution(t);
		// Make output file name for analytic solution
		const std::string name( base + std::to_string(it)+".txt");
		std::ofstream out(name);
		a.print_solution(out);
		out.close();
		t += dt;
	}
	std::cout << params.simulation_names[a.model()] << " done.\n";
	auto end = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds> (end - start);
	return duration.count();
}



#endif /* SRC_LINEAR_ANALYTIC_HH_ */
