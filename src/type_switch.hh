/*
 * type_switch.hh
 *
 *  Created on: 28 Jan 2015
 *      Author: jurak
 */

#ifndef SRC_TYPE_SWITCH_HH_
#define SRC_TYPE_SWITCH_HH_

#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/newton/newton.hh>

template <typename IGO, typename LS, typename U, bool nonlinear>
class IfNonlinear{
public:
  typedef Dune::PDELab::Newton<IGO,LS,U> PDESOLVER;
};

template <typename IGO, typename LS, typename U>
class IfNonlinear<IGO,LS,U,false>{
public:
  typedef Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,U> PDESOLVER;
};




#endif /* SRC_TYPE_SWITCH_HH_ */
