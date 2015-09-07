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

#include "parameters.hh"
#include "mgf.hh"

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

///** Boundary value function. Depends only on time.  */
//double g(double t){
//  return 0.3 + 0.2* std::sin(2*M_PI*t);
//}
///** Derivative of the boundary value function.  */
//double dg_dt(double t){
//  return  2*M_PI*0.2* std::cos(2*M_PI*t);
//}
//double dg_dt_shifted(double v, double t){ return dg_dt(t-v*v); }
/**
 * Class holding the integrand for calculation of the integral
 * (boundary layer function)
 *   \f[
 *   Z(\xi,t) = (2/\sqrt{\pi}) \int_{x/(2\delta\sqrt{t})} (g(t-\xi^2/(4 v^2)) -g(0))exp(-v^2) dv.
 *   \f]
 *   where \f$\xi\f$ is  a boundary value variable, that is, at $x=0$
 *   \f[
 *   \xi = x/\overline{\delta},\quad \overline{\delta} = \delta \sqrt{\frac{k_m\overline{\alpha}_m}{\Phi_m}}.
 *   \f]
 * For non wetting phase flux \f$Q_n\f$ we have the formula
 * \f[
 *     Q_n = \frac{4}{\sqrt{\pi}}\overline{\delta}\Phi_m \int_0^{\sqrt{t}} g'(t-\lambda^2)\, d\lambda.
 * \f]
 */
template <typename Params>
class Integrand{
  public:
    Integrand(Params const & params) : params_(params){
    	g = params.bdry_fun();
    	double delta = params.delta;
    	double perm = params.k;
    	double poro = params.poro;
    	double mean_alpha = params.mean_alpha;
    	scaled_delta = delta *std::sqrt(perm*mean_alpha/poro);
    	assert(scaled_delta > 0.0);
        // factor = 2.0/std::sqrt(M_PI);
    	factor_der = 2 * poro * scaled_delta * factor;
    	factor_1 = 2 * factor * poro * scaled_delta /params.L;
        check();

        double tend = params.tend;
        // live 10 % of space since a number of points generated is not exactly default_table_size
        dt_table_ = tend/(0.9*default_table_size_);
        last_time_index_ = 0;
        time_.resize(default_table_size_);
        tau_time_.resize(default_table_size_);
        g_tau_time_.resize(default_table_size_);
        time_(0) = 0.0;
        tau_time_(0) = 0.0;
        g_tau_time_(0) = g(0.0);
        integrate_alpha_bdry(tend);
    }
    /** Calculate the flux in the case of constant linearisation. The result is
     *  given on an equdistant time mesh and is stored in lin_flux variable.
     */
    void calculate_linear_const_flux();
    /** Print the flux calculated in calculate_linear_const_flux() on a given stream.
     *  If the flux is not already calculated, the function will call calculate_linear_const_flux().
     */
    void print_flux(std::ostream & out);
    /** Calculate the solution of constant linearization in the time instant time.
     *  The solution is stored in lin_const_solution variable.
     */
    void calculate_linear_const_solution(double time);
    /** Print the solution of constant linearization calculated by calculate_linear_const_solution()
     *  method. We do not verify that the solution is indeed calculated.
     */
    void print_solution(std::ostream & out);
    /** Calculate the solution of variable linearization in the time instant time.
     *  The solution is stored in lin_var_solution variable.
     */
    void calculate_linear_var_solution(double time);
    void calculate_linear_var_flux();
     /** Print tau(t) in file for debugging. */
    void print_tau(std::ostream & out);
    /** tau( t ). */
    double tau(double t) const;
    /** g( t(tau) ), where t() is the incerse function of tau. */
    double bdry_comp_t(double tau) const;
  private:
    const Params & params_;

    void set_x(double xx) { x_=xx; check();  xi_ = x_/scaled_delta; }

    /** Integrand in calculation of BL function Z(x,t). */
    double g_shifted(double v, double t) const {
      return factor * ( (g(t - xi_*xi_/(4*v*v)) - g(0.0)) * std::exp(-v*v) );
    }
    /**  */
    double g_t_shifted(double v, double t) const {
          return factor * ( (bdry_comp_t(t - xi_*xi_/(4*v*v)) - bdry_comp_t(0.0)) * std::exp(-v*v) );
    }
    double g_t_shifted_1(double lambda, double tau) const{
    	  return  bdry_comp_t(tau - lambda*lambda) - bdry_comp_t(0.0);
    }
    /** Integrand in calculation of the nonwetting phase flux, scalled with right factor. */
    double dg_dt_shifted(double v, double t){ return factor_der * dg_dt(t-v*v); }
    /** Integral lower bound. */
    double lower_bound(double t) const { return xi_/(2*std::sqrt(t));}
    /** Boundary condition at the time t. */
    double bdry(double t) const { return g(t); }
    /** Alpha-function. */
    double alpha(double S) const { return params_.alpha(S)/params_.mean_alpha; }
    /** Composition alpha( g(t) ). */
    double a_g(double t) const { return params_.alpha( g(t) )/params_.mean_alpha; }
    /** Calculate new entry in  new_time_field by calculating
     * \f[ \int_{t_0}^t \alpha_m( g(s)) ds,\f]
     * where \f$t_0\f$ is given by last_time_index_.
     */
    void integrate_alpha_bdry(double t);

    double scaled_delta;
    double x_= 0.0;
    double xi_ = 0.0;
    const double factor = 2.0/std::sqrt(M_PI);
    double factor_der = 0.0;
    double factor_1   = 0.0;
    const double h = 1.0E-7;
    double dt_table_ = 0.0;

    /** Linear model flux; pairs (t, flux(t)).  */
    std::vector<std::pair<double,double>> lin_flux;
    /** Linear constant solution; pairs (x, u(x,t)) at fixed time t.  */
    std::vector<std::pair<double,double>> lin_const_solution;
    /** Linear variable solution; pairs (x, u(x,t)) at fixed time t.  */
    std::vector<std::pair<double,double>> lin_var_solution;
    /** Las index of calculated time in new_time_field. */
    unsigned int last_time_index_ = 0;
    const unsigned int default_table_size_ = 10000;
    /** Pairs (t, tau(t)) where tau(t) = int_0^t alpha(g(s))ds. time_ holds t. */
    boost::numeric::ublas::vector<double> time_;
    /** Pairs (t, tau(t)) where tau(t) = int_0^t alpha(g(s))ds. tau_time_ holds tau(t).  */
    boost::numeric::ublas::vector<double> tau_time_;
    /** Values of g( tau(t) ). */
    boost::numeric::ublas::vector<double> g_tau_time_;

    void check(){
       if(x_ < 0.0  || scaled_delta <= 0.0)
         throw std::runtime_error("Integrand: x_ < 0.0 || scaled_delta <= 0.0");
    }

    /** Boundary function taken from input file. */
    std::function<double(double)> g;

    /** Derivative of the boundary value function.  */
    double dg_dt(double t){
    	return (g(t+h) - g(t-h))/(2*h);
    }

};

template <typename Params>
double Integrand<Params>::bdry_comp_t(double tau) const{
	const double TOL = 1E-10;
	double tau_max =  tau_time_[last_time_index_];
	if(tau > tau_max + TOL)
		throw std::runtime_error("Integrand<Params>::bdry_comp_t: tau > tau_max.");

	auto it = std::upper_bound(tau_time_.begin(), tau_time_.begin()+last_time_index_+1, tau);
	auto index = it - tau_time_.begin(); // this is upper bound
	if(index == 0)  return g_tau_time_[index];
	if(index > last_time_index_){
		// There is no bigger element, return largest value.
		// This can hapen only on the last element.
		assert(index == last_time_index_ + 1);
		return g_tau_time_[index-1];
//		std::cerr << tau << ", tau_max = " << tau_max << ", last_time_index = " << last_time_index_
//				 << ", index = " << index << std::endl;
//		throw std::logic_error("Internal error 542.");
	}
	const double g1 = g_tau_time_[index-1];
	const double g2 = g_tau_time_[index];
	const double t1 = tau_time_[index-1];
	const double t2 = tau_time_[index];
    double res = g1;
	if(std::abs(t2-t1) >= TOL)
           res += (tau-t1)*(g2-g1)/(t2-t1);
    return res;
}

template <typename Params>
double Integrand<Params>::tau(double t) const{
	const double TOL = 1E-10;
	double t_max =  time_[last_time_index_];
	if(t > t_max + TOL){
		std::cerr << "t= " << t << ", t_max = " << t_max << std::endl;
		throw std::runtime_error("Integrand<Params>::tau: t > t_max.");
	}
	auto it = std::upper_bound(time_.begin(), time_.begin()+last_time_index_+1, t);
	auto index = it - time_.begin(); // this is upper bound
	if(index == 0) throw std::logic_error("Internal error 542.");
	if(index > last_time_index_){
		// There is no bigger element, return largest value.
		// This can hapen only on the last element.
		assert(index == last_time_index_ + 1);
		return tau_time_[index-1];
	}
	const double g1 = tau_time_[index-1];
	const double g2 = tau_time_[index];
	const double t1 = time_[index-1];
	const double t2 = time_[index];
    double res = g1;
	if(std::abs(t2-t1) >= TOL)
           res += (t-t1)*(g2-g1)/(t2-t1);
    return res;
}


template <typename Params>
void Integrand<Params>::integrate_alpha_bdry(double t){
	unsigned int old_last_index_time = last_time_index_;
	double t_inf =  time_[last_time_index_];
	double integr=  tau_time_[last_time_index_];
	if(t_inf >= t){
		std::cerr << "Integrand<Params>::integrate_alpha_bdry(double t)  -- t_inf >= t!\n";
		// there is no need for integration -- this situation is possibly an error.
		return;
	}

	o2scl::funct11 f = std::bind(std::mem_fn<double(double)const>(&Integrand<Params>::a_g),
	    		                 this, std::placeholders::_1);
	o2scl::inte_qag_gsl<> inte_formula;
	double res=0.0, err=0.0;
	double intpart;
	double fractpart = std::modf((t - t_inf)/dt_table_, &intpart);
	if(fractpart >= 0.3) intpart++;

	for(unsigned int i=0; i<intpart; ++i){
		double time = t_inf + (i+1)*(t-t_inf)/intpart;
      	inte_formula.integ_err(f,t_inf,time,res,err);
     	last_time_index_++;
     	if(last_time_index_ >= time_.size()){
     		time_.resize(2*last_time_index_);
     		tau_time_.resize(2*last_time_index_);
     		g_tau_time_.resize(2*last_time_index_);
     	}
      	time_[last_time_index_] = time;
      	tau_time_[last_time_index_] = integr + res;
      	g_tau_time_[last_time_index_] = g( time );
	}
}

// just for debugging
template <typename Params>
void Integrand<Params>::print_tau(std::ostream & out){
	for(unsigned int i = 0; i <= last_time_index_; ++i)
		out << time_[i] << " " << tau_time_[i] << "  "
		    << g_tau_time_[i] << " " << bdry_comp_t(tau_time_[i])
			<< " " << tau(time_[i]) << "\n";
}

template <typename Params>
void Integrand<Params>::calculate_linear_const_flux(){
      double dt   = params_.dtout;
      int    Nsteps = params_.tend / dt;
      double x = 0.0, t = dt;
      lin_flux.resize(Nsteps+1);
      std::fill(lin_flux.begin(), lin_flux.end(), std::make_pair(0.0,0.0));

      o2scl::funct11 fprim = std::bind(std::mem_fn<double(double,double)>(&Integrand<Params>::dg_dt_shifted),
    		                           this, std::placeholders::_1, std::cref(t));

//      o2scl::inte_qag_gsl<> inte_formula;
      o2scl::inte_adapt_cern<> inte_formula;
      // Output flux -- one value for each time instant.
      // Time loop
      lin_flux[0].first = 0.0;
      lin_flux[0].second = 0.0;
      double res=0.0, err=0.0;
      for(int it=1; it <=Nsteps; ++it){
         // Calculate the flux by given formula
         inte_formula.integ_err(fprim,0.0,std::sqrt(t),res,err);
         lin_flux[it].first = t;
         lin_flux[it].second = res;
    //    std::cout << "err = " << err << "\n";
         t += dt;
      }
      return;
}

template <typename Params>
void Integrand<Params>::calculate_linear_var_flux(){
      double dt   = params_.dtout;
      int    Nsteps = params_.tend / dt;
      double x = 0.0, t = 0.0;
      double tau_t = 0.0;
	  const double h = 1E-6; // parameter for num diff
      lin_flux.resize(Nsteps+1);
      std::fill(lin_flux.begin(), lin_flux.end(), std::make_pair(0.0,0.0));

      o2scl::funct11 fprim = std::bind(std::mem_fn<double(double,double)const>(&Integrand<Params>::g_t_shifted_1),
    		                           this, std::placeholders::_1, std::cref(tau_t));
//      o2scl::inte_qag_gsl<> inte_formula;
      o2scl::inte_adapt_cern<> inte_formula;


      double res1=0.0, err1=0.0, res0=0.0, err0=0.0;
      // One sided diference for t=0
      inte_formula.integ_err(fprim,0.0,std::sqrt(0.0 + h),res1,err1);
      lin_flux[0].first = t;
      lin_flux[0].second = factor_1 * a_g(0.0) * (res1 - 0.0)/h;
      // Output flux -- one value for each time instant.
      // Time loop
      t = dt;
      for(int it=1; it <Nsteps; ++it){
         // Calculate the flux by given formula
         tau_t = tau(t)+ h/2;
         inte_formula.integ_err(fprim,0.0,std::sqrt(tau_t + h/2),res1,err1);
         tau_t = tau(t)- h/2;
         inte_formula.integ_err(fprim,0.0,std::sqrt(tau_t - h/2),res0,err0);
         lin_flux[it].first = t;
         lin_flux[it].second = factor_1 * a_g(t) * (res1 - res0)/h;
    //     std::cout << factor_1 << " " << a_g(t) << " " << (res1 - res0)/h << std::endl;
    //    std::cout << "err = " << err << "\n";
         t += dt;
      }
      tau_t = tau(t);
      inte_formula.integ_err(fprim,0.0,std::sqrt(tau_t),res1,err1);
      tau_t = tau(t)-h;
      inte_formula.integ_err(fprim,0.0,std::sqrt(tau_t - h),res0,err0);
      lin_flux[Nsteps].first = t;
      lin_flux[Nsteps].second = factor_1 * a_g(t) * (res1 - res0)/h;

      return;
}

template <typename Params>
void Integrand<Params>::print_flux(std::ostream & out){
	// Note. Flux must be calculated firs.
    for(unsigned int i = 0; i < lin_flux.size(); ++i)
    	out << lin_flux[i].first << "  " << lin_flux[i].second <<"\n";
}


template<typename Params>
void Integrand<Params>::calculate_linear_const_solution(double time) {
	double L = params_.L;
	MGF<Params> d1(params_);
	std::vector<double> x2;
	d1.double_side_interval(x2);
	lin_const_solution.resize(params_.N);

	o2scl::funct11 f = std::bind(
			std::mem_fn<double(double, double)const>(&Integrand<Params>::g_shifted),
			this, std::placeholders::_1, std::cref(time));
//	  o2scl::inte_qag_gsl<> inte_formula;
	o2scl::inte_adapt_cern<> inte_formula;

	double lb, ub, res1, res2, err;
	if (time == 0.0)
		for (int i = 0; i < params_.N; ++i) {
			lin_const_solution[i].first = x2[i];
			lin_const_solution[i].second = bdry(0.0);
		}
	else
		for (int i = 0; i < params_.N; ++i) {
			set_x(x2[i]);
			lb = lower_bound(time);  // depends on x!
			ub = std::max(6.0, lb + 1);
			// calculate Z(x,t)
			inte_formula.integ_err(f, lb, ub, res1, err);
//	        std::cout << "err = " << err <<"\n";
			set_x(L - x2[i]);
			lb = lower_bound(time);  // depends on x!
			ub = std::max(6.0, lb + 1);
			inte_formula.integ_err(f, lb, ub, res2, err);

			lin_const_solution[i].first = x2[i];
			lin_const_solution[i].second = bdry(0.0) + res1 + res2;
		}
	return;
}

template<typename Params>
void Integrand<Params>::calculate_linear_var_solution(double time) {
	double tau_t = tau(time);
	double L = params_.L;
	MGF<Params> d1(params_);
	std::vector<double> x2;
	d1.double_side_interval(x2);
	lin_var_solution.resize(params_.N);

	o2scl::funct11 f = std::bind(
			std::mem_fn<double(double, double)const>(&Integrand<Params>::g_t_shifted),
			this, std::placeholders::_1, std::cref(tau_t));
//	  o2scl::inte_qag_gsl<> inte_formula;
	o2scl::inte_adapt_cern<> inte_formula;

	double lb, ub, res1, res2, err;
	if (time == 0.0)
		for (int i = 0; i < params_.N; ++i) {
			lin_var_solution[i].first = x2[i];
			lin_var_solution[i].second = bdry(0.0);
		}
	else
		for (int i = 0; i < params_.N; ++i) {
			set_x(x2[i]);
			lb = lower_bound(tau_t);  // depends on x!
			ub = std::max(6.0, lb + 1);
			// calculate Z(x,t)
			inte_formula.integ_err(f, lb, ub, res1, err);
//	        std::cout << "err = " << err <<"\n";
			set_x(L - x2[i]);
			lb = lower_bound(tau_t);  // depends on x!
			ub = std::max(6.0, lb + 1);
			inte_formula.integ_err(f, lb, ub, res2, err);

			lin_var_solution[i].first = x2[i];
			lin_var_solution[i].second = bdry(0.0) + res1 + res2;
		}
	return;
}


template <typename Params>
void Integrand<Params>::print_solution(std::ostream & out){
	if(params_.model == Params::analytic_const)
	for(unsigned int i = 0; i < lin_const_solution.size(); ++i)
	    	out << lin_const_solution[i].first << "  " << lin_const_solution[i].second <<"\n";
	else
	for(unsigned int i = 0; i < lin_var_solution.size(); ++i)
	    	out << lin_var_solution[i].first << "  " << lin_var_solution[i].second <<"\n";

}

// compilation
//  g++ -std=c++11 -g linear_analytic.cc  -o analytic /usr/local/lib/libo2scl.a -lgsl -lblas
template<typename Params>
void lin_analytic(Params params) {
	Integrand<Params> a(params);

	double dt = params.dtout;
	int Nsteps = params.tend / dt;
	double t = 0.0;

	const std::string base( params.str_sname + params.simulation_names[params.model] );
    const std::string tau = base + "-tau.txt";
    std::ofstream out_tau(tau);
    a.print_tau(out_tau);
    out_tau.close();

    if(params.model == Params::analytic_const)
	     a.calculate_linear_const_flux();
    else
	     a.calculate_linear_var_flux();
	// Output flux -- one value for each time instant.
	const std::string flux( base + "-flux.txt");
	std::ofstream out_flux(flux);
	a.print_flux(out_flux);
	out_flux.close();
	// Time loop
	for (int it = 0; it <= Nsteps; ++it) {
		if(params.model == Params::analytic_const)
		     a.calculate_linear_const_solution(t);
		else
		     a.calculate_linear_var_solution(t);
		// Make output file name for analytic solution
		const std::string name( base + std::to_string(it)+".txt");
		std::ofstream out(name);
		a.print_solution(out);
		out.close();
		t += dt;
	}
	return;
}



#endif /* SRC_LINEAR_ANALYTIC_HH_ */
