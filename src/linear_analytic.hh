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

#include <o2scl/inte_qag_gsl.h>

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
    Integrand(Params const & params){
    	g = params.bdry_fun();
    	double delta = params.delta;
    	double perm = params.k;
    	double poro = params.poro;
    	double mean_alpha = params.mean_alpha;
    	scaled_delta = delta *std::sqrt(perm*mean_alpha/poro);
    	assert(scaled_delta > 0.0);
    	factor_der = 2 * poro * scaled_delta * factor;
        check();
    }
    void set_x(double xx) { x_=xx; check();  xi_ = x_/scaled_delta; }

    /** Integrand in calculation of BL function Z(x,t). */
    double g_shifted(double v, double t){
      return factor * ( (g(t - xi_*xi_/(4*v*v)) - g(0.0)) * std::exp(-v*v) );
    }
    /** Integrand in calculation of the nonwetting phase flux, scalled with right factor. */
    double dg_dt_shifted(double v, double t){ return factor_der * dg_dt(t-v*v); }
    /** Integral lower bound. */
    double lower_bound(double t) const { return xi_/(2*std::sqrt(t));}
    double xi() const { return xi_; }
    /** x coordinate. */
    double get_x() const { return x_; }
    /** Time t. */
    double bdry(double t) const { return g(t); }
    double scaled_delta;
  private:
    double x_= 0.0;
    double xi_ = 0.0;
    const double factor = 2.0/std::sqrt(M_PI);
    double factor_der = 0.0;
    const double h = 1.0E-7;
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


// compilation
//  g++ -std=c++11 -g linear_analytic.cc  -o analytic /usr/local/lib/libo2scl.a -lgsl -lblas
template <typename Params>
void lin_analytic(Params params){
//   double M = find_upper_bound(1.0E-6);
//   std::cout << "TOL = 1E-6, M = " << M << std::endl;
//   M = find_upper_bound(1.0E-8);
//   std::cout << "TOL = 1E-8, M = " << M << std::endl;
//   M = find_upper_bound(1.0E-10);
//   std::cout << "TOL = 1E-10, M = " << M << std::endl;
//   M = find_upper_bound(1.0E-12);
//   std::cout << "TOL = 1E-12, M = " << M << std::endl;
//   M = find_upper_bound(1.0E-14);
//   std::cout << "TOL = 1E-14, M = " << M << std::endl;
//
  double L    = params.L;
  double dt   = params.dtout;
  int    Nsteps = params.tend / dt;

  MGF<Params>  d1(params);
  std::vector<double> x2;
  d1.double_side_interval(x2);

  double x = 0.0, t = dt;
  Integrand<Params> a(params);

  o2scl::funct11 f     = std::bind(std::mem_fn<double(double,double)>(&Integrand<Params>::g_shifted),
		                           &a, std::placeholders::_1, std::cref(t));
  o2scl::funct11 fprim = std::bind(std::mem_fn<double(double,double)>(&Integrand<Params>::dg_dt_shifted),
		                           &a, std::placeholders::_1, std::cref(t));

  o2scl::inte_qag_gsl<> inte_formula;
  // Output flux -- one value for each time instant.
  std::string flux(params.str_sname+params.simulation_names[params.model]+"-flux.txt");
  std::ofstream out_flux(flux);
  // Time loop
  for(int it=1; it <=Nsteps; ++it){
     // Calculate the flux by given formula
     double res=0.0, err=0.0;
     inte_formula.integ_err(fprim,0.0,std::sqrt(t),res,err);
     out_flux << t << " " << res << "\n";
//    std::cout << "err = " << err << "\n";


     // Make output file name for analytic solution
     std::string name = params.str_sname + params.simulation_names[params.model];
   //  std::string name(params.date_and_time+"/anLin-");
     name +=std::to_string(it);
     name += ".txt";
     std::ofstream out(name);
     // Here we calculate analytic solution in the whole domain (0,L).
     // However, we calculate only the boundary layer
     for(int i=0; i<params.N; ++i)
     {
        a.set_x(x2[i]);
        double lb = a.lower_bound(t), ub = std::max(6.0, lb+1), res1, res2, err;
        // calculate Z(x,t)
        inte_formula.integ_err(f,lb,ub,res1,err);
        a.set_x(L - x2[i]);
        lb = a.lower_bound(t);
        ub = std::max(6.0, lb+1);
        inte_formula.integ_err(f,lb,ub,res2,err);
        out << a.get_x() << "   " << a.bdry(0.0)+res1+res2 <<"\n";
     }
     out.close();
     t += dt;
  }
  out_flux.close();
  return;
}



#endif /* SRC_LINEAR_ANALYTIC_HH_ */
