/*
 * test.hh
 *
 *  Created on: 14 Jan 2015
 *      Author: jurak
 */

#ifndef SRC_TEST_HH_
#define SRC_TEST_HH_

// grid function from analytic function. We need to provide
// only evaluateGlobal and the base class will provide
// evaluate(element, localCoo, value)

template <typename A, typename B>
using AGFB = Dune::PDELab::AnalyticGridFunctionBase<A,B>;

template <typename GV, typename RF, int N>
using AGFT = Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,N>;

template<typename GV, typename RF>
class Square : public AGFB< AGFT<GV,RF,1>, Square<GV,RF> > {
public:
  typedef AGFT<GV,RF,1> Traits;
  typedef AGFB<Traits,Square<GV,RF> > B;

  Square(const GV& gv) : B(gv) {}
  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {
//    const double pi = M_PI;
    y = x*x;
  }
};




#endif /* SRC_TEST_HH_ */
