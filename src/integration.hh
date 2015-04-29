/*
 * integration.hh
 *
 *  Created on: 9 Jan 2015
 *      Author: jurak
 */

#ifndef SRC_INTEGRATION_HH_
#define SRC_INTEGRATION_HH_

#include <dune/geometry/quadraturerules.hh>



/**
 * Routine for volume integration of a grid function
 *  DGF = DuneGridFunction
 *  @param dgf = DuneGridFunction to integrate
 *  @param order = order of the integration formula
 */
template<typename DGF>
  double
  volume_integral (DGF const & dgf, int order)
  {
    double integral = 0.0;
    auto const & gv = dgf.getGridView ();
    const int dim = gv.dimension;

    // Iterate by all elements
    auto el_it = gv.template begin<0> ();
    for (; el_it != gv.template end<0> (); ++el_it)
      {
	const auto & gtype = el_it->geometry ().type ();
	const auto & rule = Dune::QuadratureRules<double, dim>::rule (gtype, order);

	double elemIntegral = 0.0;
	for (auto ii = rule.begin (); ii != rule.end (); ++ii)
	  {
	    const auto & xi = ii->position ();
	    double omegai = ii->weight ();
	    double detJaci = el_it->geometry ().integrationElement (xi);
	    Dune::FieldVector<double, 1> functioni = 0.0;
	    dgf.evaluate(*el_it, xi, functioni);
	    elemIntegral += functioni * omegai * detJaci;
	  }
	integral += elemIntegral;
      }
    return integral;
  }

template<typename DGF>
  double
  boundary_integral (DGF const & dgf, int order)
  {
    double integral = 0.0;
    auto const & gv = dgf.getGridView ();
    const int dim = gv.dimension;

    // By all elements
    auto el_it = gv.template begin<0> ();
    for (; el_it != gv.template end<0> (); ++el_it)
      {
	const auto el_geo = el_it->geometry();
	auto isit_end = gv.iend(*el_it);
	auto isit     = gv.ibegin(*el_it);

	   // By all sides
	   for( ; isit != isit_end; ++isit){
	       if(isit->boundary()){
		   auto outerNormal = isit->centerUnitOuterNormal();
		   const auto igeo = isit->geometry();  // intersection geometry
		   const auto gt = igeo.type();
		   const auto& rule = Dune::QuadratureRules<double, dim-1>::rule(gt, order);

		   double side_integral = 0.0;

		   auto iq = rule.begin();
		   for (; iq != rule.end(); ++iq) {
		       typename DGF::Traits::RangeType fval;
		       dgf.evaluate(*el_it, el_geo.local( igeo.global(iq->position())), fval );
		       double weight = iq->weight();
		                      // | det (grad g) |
		       double detjac = igeo.integrationElement(iq->position());
		       side_integral += (fval * outerNormal) * weight * detjac;
		   }
		   integral += side_integral;
	       }
	   }// end by all sides
      }// end by all elements
    return integral;
  }



#endif /* SRC_INTEGRATION_HH_ */
