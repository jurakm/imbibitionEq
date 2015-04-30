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

// Class that calculates integral of given function in a
// form of a table, ad offers function that interpolates the table,
// Prerequisite:
//      Given function is a member function of class FunctionClass.
//      It s defined on [0,1] with the range in R.
template <typename FunctionClass>
class Table{
public:
  using FunctionClassPtr = double (FunctionClass::*)(double) const;

  Table() : N(0), h(0.0), pobject(nullptr), pfun(nullptr){}
  Table(FunctionClass const * pobject_, FunctionClassPtr pfun_, unsigned int N_) :
    N(0), h(0.0), pobject(nullptr), pfun(nullptr)
  { init(pobject_, pfun_, N_);}
  double interpolate(double sw) const;
  unsigned int NofPts() const { return N; }
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

// Purely artificial alpha() with one parameter.
class ArtifImbibitionFunctions{
public:
  ArtifImbibitionFunctions() : a(0.0){}
  ArtifImbibitionFunctions(double a_, unsigned int N_=Ndefault) : a(a_),
  		       table(this, &ArtifImbibitionFunctions::alpha, N_) {}
  void init(double a_, unsigned int N_=Ndefault)
  {
    a = a_;
    table.init(this, &ArtifImbibitionFunctions::alpha,N_);
  }

  double alpha(double swe) const { return a*swe*(1-swe); }
  double beta(double swe) const { return table.interpolate(swe); }
  unsigned int NofPts() const { return table.NofPts(); }
private:
  double a;
  Table<ArtifImbibitionFunctions> table;
  static const unsigned int Ndefault = 10000;
};



template <typename TwoPhaseFunctions>
class RealImbibitionFunctions{
public:
  using Params = typename TwoPhaseFunctions::Params;
  RealImbibitionFunctions(){}
  RealImbibitionFunctions(Params const & params_, double muw_, double mun_, unsigned int N_= Ndefault)
  {
    init(params_, muw_, mun_, N_);
  }
  void init(Params const & params_, double muw_, double mun_, unsigned int N_= Ndefault){
    params = params_; muw = muw_; mun = mun_;
    table.init(this, &RealImbibitionFunctions<TwoPhaseFunctions>::alpha, N_);
  }
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
  double beta(double swe) const { return table.interpolate(swe); }
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
  beta_values.resize(N+1);
  for(auto & x : beta_values) x = 0.0;
  beta_values[0] = 0.0;
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
  double fraction = std::modf(sw/h, &int_part);
  unsigned int i = int_part;

  double val = fraction*beta_values[i+1] + (1-fraction)*beta_values[i];
  return val;
}

// ----------------------------  imbibitionFunctions --------------------


template <typename ImbibitionFunctions>
void plot(std::string const & file_name, ImbibitionFunctions const & obj, unsigned int n)
{
    if(n == 0) throw std::runtime_error("n=0");
    double hh = 1.0/n;

    std::ofstream file (file_name);
    file << "    x         alpha(x)            beta(x)\n";
    for (unsigned int i = 0; i <= n; ++i)
      {
	const double xi = i*hh;
	file << std::setw(10) << std::setprecision(8) << xi << "   "
	     << std::setw(16) << std::setprecision(12) << obj.alpha(xi)  << "   "
	     << std::setw(16) << std::setprecision(12) << obj.beta(xi) << "\n";
      }
    file.close ();
}

#endif
