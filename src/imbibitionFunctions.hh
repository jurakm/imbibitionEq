#ifndef __IMBIBITION_FUNCTIONS_IS_INCLUDED__
#define __IMBIBITION_FUNCTIONS_IS_INCLUDED__
/*
 * imbibitionFunctions.hh
 *
 *  Created on: 24 Apr 2015
 *      Author: jurak
 */
#include <vector>
#include <stdexcept>
#include <string>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "tableXY.hh"


// TODO Put this class in use -- replace Table<>
template <typename FunctionClass>
class TabularFunction : private TableXY<std::vector<double> >{
public:
  using FunctionClassPtr = double (FunctionClass::*)(double) const;

  /// Default constructor constructs an empty table.
  TabularFunction() : TableXY<std::vector<double> >(), pobject(nullptr), pfun(nullptr){}

 /// Constructor. Actual initialization of the table is done in function init().
 /// @param pobject_ FunctionClass object.
 /// @param pfun_    Pointer to a FunctionClass method.
 /// @param N_       Number of points in the table.
  TabularFunction(FunctionClass const * pobject_, FunctionClassPtr pfun_, unsigned int N_=10000) :
	  TableXY<std::vector<double> >(N_), pobject(nullptr), pfun(nullptr)
  { init(pobject_, pfun_, N_);}
  /// Return the function value by linear interpolation.
  /// @param sw saturation value in [0,1].
  double interpolate(double sw) const;
  double interpolate_inverse(double u) const;
  /// Number of points in the table
  unsigned int NofPts() const { return last()+1; }
  /// Construction of the function table. If the table is already in use it will be cleared first.
 /// @param pobject_ FunctionClass object.
 /// @param pfun_    Pointer to a FunctionClass method.
 /// @param N_       Number of points in the table.
  void init(FunctionClass const * pobject_, FunctionClassPtr pfun_, unsigned int N_);
private:

  FunctionClass const * pobject;
  double (FunctionClass::*pfun)(double) const;
};


template <typename FunctionClass>
inline void TabularFunction<FunctionClass>::init(FunctionClass const * pobject_,
				       FunctionClassPtr pfun_, unsigned int N_){
  resize(N_ + 1); // N_ intervals means N_+1 points.
  double h = 0.0;
  pobject = pobject_;
  pfun = pfun_;

  if(N_ == 0) throw std::runtime_error("N=0");
  h = 1.0/N_;  // Domain is [0,1].  TODO make domain more general.

  // If the table was already in use this will clear it.
  clean();
  push_back(0.0,0.0); // first element in the table

//  beta_values[0] = 0.0;
  for(unsigned int i = 0; i < N_; ++i){
      // calculate integral on (x_i,x_{i+1})
      const double xm = i*h;
      const double xp = (i+1)*h;
      const double xc = (xm+xp)/2;
      // Simpson's rule
      const double integral = (h/6)*( (pobject->*pfun)(xm) + 4*(pobject->*pfun)(xc) + (pobject->*pfun)(xp) );
//      std::cout << (pobject->*pfun)(xm) << std::endl;
      const double prev_int = get_y(i);
      push_back(xp, prev_int + integral);
  }
}

template <typename FunctionClass>
double TabularFunction<FunctionClass>::interpolate(double sw) const
{
  // Extension by constant outside [0,1].
  if(sw < 0.0)  sw = 0.0;
  if(sw > 1.0)  sw = 1.0;

  const double val = interpolate_y(sw);
  return val;
//  double int_part = 0.0;
//  const double fraction = std::modf(sw/h, &int_part);
//  const unsigned int i = static_cast<unsigned int>(int_part);
//
//  const double val = fraction * beta_values[i+1] + (1-fraction) * beta_values[i];
//  return val;
}

template <typename FunctionClass>
double TabularFunction<FunctionClass>::interpolate_inverse(double u) const
{
  // Extension by constant outside [0,max].
  if(u < 0.0)  u = 0.0;
  if(u > get_last_y())  u = get_last_y();

  const double val = interpolate_x(u);
  return val;
}


/// @brief Class that calculates integral of a given function in a form of a table,
///         and offers function that interpolates the table. Function that is to integrated
///         is given as member function of the class given as the template parameter.
/// Prerequisite:
///      Given function is a member function of class FunctionClass.
///      It s defined on [0,1] with the range in R.
 /// @tparam FunctionClass Class of functions.
template <typename FunctionClass>
class Table{
public:
  using FunctionClassPtr = double (FunctionClass::*)(double) const;

  /// Default constructor constructs an empty table.
  Table() : N(0), h(0.0), pobject(nullptr), pfun(nullptr){}

 /// Constructor. Actual initialization of the table is done in function init().
 /// @param pobject_ FunctionClass object.
 /// @param pfun_    Pointer to a FunctionClass method.
 /// @param N_       Number of points in the table.
  Table(FunctionClass const * pobject_, FunctionClassPtr pfun_, unsigned int N_) :
    N(0), h(0.0), pobject(nullptr), pfun(nullptr)
  { init(pobject_, pfun_, N_);}
  /// Return the function value by linear interpolation.
  /// @param sw saturation value in [0,1].
  double interpolate(double sw) const;
  /// Number of points in the table
  unsigned int NofPts() const { return N; }
  /// Construction of the function table. If the table is already in use it will be cleared first.
 /// @param pobject_ FunctionClass object.
 /// @param pfun_    Pointer to a FunctionClass method.
 /// @param N_       Number of points in the table.
  void init(FunctionClass const * pobject_, FunctionClassPtr pfun_, unsigned int N_);
private:
  unsigned int N;
  double h;

  FunctionClass const * pobject;
  double (FunctionClass::*pfun)(double) const;
  // Table of values for integral of the function *pfun.
  std::vector<double> beta_values;
  // Make a table of values
};

// Virtual base for the functions alpha and beta.
class ImbibitionFunctions{
public:
	virtual double alpha(double swe) const = 0;
	virtual double beta (double swe) const = 0;
//	virtual double gamma(double u)   const = 0;
	virtual double pc   (double swe) const = 0;
	virtual double sw   (double pc ) const = 0;
	virtual unsigned int NofPts() const = 0;

	virtual ~ImbibitionFunctions() {}
};
/// \brief Class giving nonlinear diffusion coefficient and its integral.
///
/// Purely artificial \f$ \alpha(S)= aS(1-S) \f$ function with one parameter.
/// Function \f$\beta(S) = \int_0^S \alpha(u) du\f$ is given in a form of a table
/// formed by numerical integration (composite Simpson's rule).
class ArtifImbibitionFunctions : public ImbibitionFunctions{
public:
	/// Default constructor. init() must be called before any use.
  ArtifImbibitionFunctions() : a(0.0){}
    /// Constructor.
	/// @param a_ = coefficient a in \f$ \alpha(S)= aS(1-S) \f$
    /// @param N_ = number of points in the table for \f$\beta(S) = \int_0^S \alpha(u) du.\f$
  ArtifImbibitionFunctions(double a_, unsigned int N_=Ndefault) : a(a_),
  		       table(this, &ArtifImbibitionFunctions::alpha, N_) {}
  void init(double a_, unsigned int N_=Ndefault)
  {
    a = a_;
    table.init(this, &ArtifImbibitionFunctions::alpha,N_);
  }
  /// \f$ \alpha(S)= aS(1-S) \f$
  double alpha(double swe) const { return a*swe*(1-swe); }
  ///  \f$\beta(S) = \int_0^S \alpha(u) du.\f$
  double beta(double swe) const { return table.interpolate(swe); }
  /// dummy capillary pressure
  double pc   (double swe) const {return swe; }
  /// dummy inverse of the capillary pressure
  double sw   (double pc ) const {return pc; }
  /// Number of points in the table table for \f$\beta(S) = \int_0^S \alpha(u) du.\f$
  unsigned int NofPts() const { return table.NofPts(); }
private:
  double a;
  Table<ArtifImbibitionFunctions> table;
  static const unsigned int Ndefault = 10000;
};


/// \brief Class giving nonlinear diffusion coefficient and its integral.
///
/// Function  \f$ \alpha(S)= - \lambda_w(S) \lambda_n(S) p_c'(S)/\lambda(S) \f$
/// based on two-phase flow  functions.
/// Function \f$\beta(S) = \int_0^S \alpha(u) du\f$ is given in a form of a table
/// formed by numerical integration (composite Simpson's rule).
/// @tparam = Two-phase flow function class.
template <typename TwoPhaseFunctions>
class RealImbibitionFunctions : public ImbibitionFunctions{
public:
  using Params = typename TwoPhaseFunctions::Params;
  /// Default constructor. init() must be called before the use.
  RealImbibitionFunctions(){}
  /// Constructor.
  /// @param params_ = parameters for two-phase flow functions which are given through template
  /// parameter.
  /// @param muw_ = viscosity of the wetting phase.
  /// @param mun_ = viscosity of the non wetting phase.
  /// @param N_ = number of points in the table.
  RealImbibitionFunctions(Params const & params_, double muw_, double mun_, unsigned int N_= Ndefault)
  {
    init(params_, muw_, mun_, N_);
  }
  /// Actuall constructor.
  /// @param params_ = parameters for two-phase flow functions which are given through template
  /// parameter.
  /// @param muw_ = viscosity of the wetting phase.
  /// @param mun_ = viscosity of the non wetting phase.
  /// @param N_ = number of points in the table.
  void init(Params const & params_, double muw_, double mun_, unsigned int N_= Ndefault){
    params = params_; muw = muw_; mun = mun_;
    table.init(this, &RealImbibitionFunctions<TwoPhaseFunctions>::alpha, N_);
  }

  /// \f$ \alpha(S)\f$
  double alpha(double swe) const{
    // extension by zero outside of the [0,1]
    if(swe <= 0.0 or swe >= 1.0) return 0.0;
     const double lambdaw = TwoPhaseFunctions::krw(params, swe) / muw;
     const double lambdan = TwoPhaseFunctions::krn(params, swe) / mun;
     const double lambda = lambdaw + lambdan;
     const double a =  lambdaw * lambdan / lambda;
     const double adpc = -a * TwoPhaseFunctions::dpc_dsw(params, swe);
     return adpc;
  }
  ///  \f$\beta(S) = \int_0^S \alpha(u) du.\f$
  double beta(double swe) const { return table.interpolate(swe); }
  /// the capillary pressure
  double pc   (double swe) const {return TwoPhaseFunctions::pc(params, swe); }
  /// inverse of the capillary pressure
  double sw   (double pc ) const {return TwoPhaseFunctions::sw(params, pc); }
  /// Number of points in the table table for \f$\beta(S) = \int_0^S \alpha(u) du.\f$
  unsigned int NofPts() const { return table.NofPts(); }
private:
  Params params;
  double muw = 0.0, mun = 0.0; // viscosity of wetting and non-wetting phase
  static const unsigned int Ndefault = 10000;

  Table< RealImbibitionFunctions<TwoPhaseFunctions> > table;
};


// -----------------  Table -----------------------------------

template <typename FunctionClass>
inline void Table<FunctionClass>::init(FunctionClass const * pobject_,
				       FunctionClassPtr pfun_,
				       unsigned int N_){
  N = N_;
  h = 0.0;
  pobject = pobject_;
  pfun = pfun_;

  if(N == 0) throw std::runtime_error("N=0");
  h = 1.0/N;
//  std::cout << "init called with N = " << N << std::endl;
  // If the table was already in use this will clear it.
  beta_values.resize(N+1);
  for(auto & x : beta_values) x = 0.0;

//  beta_values[0] = 0.0;
  for(unsigned int i = 0; i < N; ++i){
      // calculate integral on (x_i,x_{i+1})
      const double xm = i*h;
      const double xp = (i+1)*h;
      const double xc = (xm+xp)/2;
      // Simpson's rule
      const double integral = (h/6)*( (pobject->*pfun)(xm) + 4*(pobject->*pfun)(xc) + (pobject->*pfun)(xp) );
//      std::cout << (pobject->*pfun)(xm) << std::endl;
      beta_values[i+1] = beta_values[i] + integral;

  }
}

template <typename FunctionClass>
double Table<FunctionClass>::interpolate(double sw) const
{
  // Extension by constant outside [0,1].
  if(sw <= 0.0)  sw = 0.0;
  if(sw >= 1.0)  sw = 1.0;

  double int_part = 0.0;
  const double fraction = std::modf(sw/h, &int_part);
  const unsigned int i = static_cast<unsigned int>(int_part);

  const double val = fraction * beta_values[i+1] + (1-fraction) * beta_values[i];
  return val;
}



#endif
