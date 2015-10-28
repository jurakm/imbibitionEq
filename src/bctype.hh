#ifndef __BCTYPE_HH_IS_INCLUDED__
#define __BCTYPE_HH_IS_INCLUDED__

#include <cmath>
#include "parameters.hh"


/** \brief constraint parameter class selecting boundary condition type */
class BCTypeParam
  : public Dune::PDELab::DirichletConstraintsParameters 
{
  double time;
public:
  template<typename I>
  bool isDirichlet(const I & intersection,   
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    return true;
  }

  //! set time for subsequent evaluation
  void setTime (double t) { time = t; }
  double getTime() const {return time; }
};

/** \brief  The Dirichlet boundary condition.
 *   It gives also the initial condition at  t=0.
 */
template<typename GV, typename RF, typename Params>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, BCExtension<GV,RF,Params> >
{
  const GV& gv;
  RF time;
  Params const & params;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  BCExtension (const GV& gv_, Params const & params_) : gv(gv_), params(params_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
//       const double TOL = 1e-10;
//       auto xglobal = e.geometry().global(xlocal);
//       double L = params.L;
//       int dim = GV::dimension;
//       if(time == 0.0) y = 0.0; // initial condition
//       else            y = params.bdry(time);
// 
//       y = params.bdry(0.0);  // hack
// 
//       for(int d = 0; d < dim; ++d){
//            if(xglobal[d] < TOL or xglobal[d] > L - TOL) // this is the boundary condition
//               y = params.bdry(time);
//       }
       y = params.bdry(time);
      return;
  }

  inline const GV& getGridView () {return gv;}

  // Set time before calling evaluate().
  void setTime (double t) {time = t;}
};
#endif
