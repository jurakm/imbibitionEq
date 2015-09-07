/*
 * integration.hh
 *
 *  Created on: 9 Jan 2015
 *      Author: jurak
 */

#ifndef SRC_INTEGRATION_HH_
#define SRC_INTEGRATION_HH_

#include <vector>

#include <dune/geometry/quadraturerules.hh>

/**  @brief Class for calculation of volume and boundary integrals.
 *
 *
 *
 */
class Integration{
public:
	/// Constructor
	Integration(int vol_order_ = 3, int bdr_order_ = 3) : vol_order(vol_order_), bdr_order(bdr_order_)
    {
		volume_values.reserve(1024); bdry_values.reserve(1024); // for efficiency
    }
	/// Calculates volume and boundary integrals.
	/// It calls private methods volume_integral() and boundary_integral()
	template <typename DGF, typename DGFGrad>
	void integrate(double time, DGF const & udgf, DGFGrad const & grad_udgf);
	/// Calculates the time derivative of the volume integral.
	void volume_derivative();
	/// Prints the results of integration to a file.
    template <typename Params>
	void print(std::string const & file_name, Params const & params);
private:
	/// Calculates volume integral.
	template<typename DGF>
	double volume_integral(DGF const & dgf, int order);

	/// Calculates boundary integral.
	template<typename DGF>
	double boundary_integral(DGF const & dgf, int order);

	std::vector<std::pair<double, double> > volume_values; // (t, int (t))
	std::vector<std::pair<double, double> > bdry_values;  // (t, int (t))
	int vol_order; ///< Integration order in volume integration
	int bdr_order; ///< Integration order in boundary integration
	std::vector<std::pair<double, double> > volume_values_der;
};


template<typename Params>
void Integration::print(std::string const & file_name, Params const & params) {
	std::ofstream out1(file_name);

	out1 << "#     t        Phi d/dt int S   k delta^2 alpha(g(t)) int bdry grad S .n    bdry(t)\n";
	const double kdd = params.k * params.delta * params.delta;
//	std::cout << "kdd = " << kdd << std::endl;
//	std::cout << "params.mean_alpha = " << params.mean_alpha << std::endl;
	for (unsigned int i = 0; i < volume_values.size(); ++i) {
		double time = volume_values_der[i].first;
		if (std::abs(time - bdry_values[i].first) >= 1.0E-12)
			throw std::runtime_error(
					std::string("Time error! i = ") + std::to_string(i) + " "
							+ std::to_string(time) + " "
							+ std::to_string(bdry_values[i].first));

		out1 << std::setw(10) << std::setprecision(4) << time
				<< " "  // time
				<< std::setw(12) << std::setprecision(6)
				<< params.poro * volume_values_der[i].second; //  " = d/dt int S"
		if (params.model == Params::nonlinear)
			out1 << "                 " << std::setw(12) << std::setprecision(6)
					<< params.alpha(params.bdry(time)) * kdd
							* bdry_values[i].second;
		// "= k delta^2 alpha(g(t)) int bdry grad S .n  "
		else if (params.model == Params::new_nonlinear)
					out1 << "                 " << std::setw(12) << std::setprecision(6)
							<<  kdd * bdry_values[i].second;
				// "= k delta^2 int bdry grad beta(S) .n  "
		else if (params.model == Params::variable_linear)
			out1 << "                 " << std::setw(12) << std::setprecision(6)
					<< params.alpha_reg(params.bdry(time)) * kdd * bdry_values[i].second;
		// "= k delta^2 alpha(g(t)) int bdry grad S .n  "
		else if (params.model == Params::constant_linear)
			out1 << "                 " << std::setw(12) << std::setprecision(6)
					<< params.mean_alpha * kdd * bdry_values[i].second;
		// "= k delta^2 mean_alpha int bdry grad S .n  "

		out1 << "             " << std::setw(12) << std::setprecision(6)
				<< params.bdry(time) << "\n";
	}
	out1.close();
}

void Integration::volume_derivative(){
	volume_values_der.resize(volume_values.size()); // (t, d/dt int S(t))

		double dt0 = volume_values[1].first - volume_values[0].first;
		double dS0 = volume_values[1].second - volume_values[0].second;
		volume_values_der[0] = std::make_pair(volume_values[0].first, dS0 / dt0);

		unsigned int nn = volume_values.size() - 1;
		for (unsigned int i = 1; i < nn; ++i) {
			double dt = volume_values[i + 1].first - volume_values[i - 1].first;
			double dS = volume_values[i + 1].second - volume_values[i - 1].second;

			volume_values_der[i] = std::make_pair(volume_values[i].first, dS / dt);
		}
		double dtn = volume_values[nn].first - volume_values[nn - 1].first;
		double dSn = volume_values[nn].second - volume_values[nn - 1].second;
		volume_values_der[nn] = std::make_pair(volume_values[nn].first, dSn / dtn);
}

template <typename DGF, typename DGFGrad>
void Integration::integrate(double time, DGF const & udgf, DGFGrad const & grad_udgf){
	// volume integral of saturation
	double intOfsol = volume_integral(udgf, vol_order);
	volume_values.push_back(std::make_pair(time, intOfsol));
//	      std::cout << "(t, int(u)) = (" << time <<","
//	                 << std::setprecision(12) << intOfsol << ")\n";

	// boundary integral of the gradient of saturation/beta(S)
	double intOverBoundary = boundary_integral(grad_udgf, bdr_order);
	bdry_values.push_back(std::make_pair(time, intOverBoundary));
	//      std::cout << "(t, int_bdry(grad u)) = (" << time <<","
	//                << std::setprecision(12) << intOverBoundary << ")\n";
}



/**
 * Routine for volume integration of a grid function
 *  DGF = DuneGridFunction
 *  @param dgf = DuneGridFunction to integrate
 *  @param order = order of the integration formula
 */
template<typename DGF>
double Integration::volume_integral(DGF const & dgf, int order) {
	double integral = 0.0;
	auto const & gv = dgf.getGridView();
	const int dim = gv.dimension;

	// Iterate by all elements
	auto el_it = gv.template begin<0>();
	for (; el_it != gv.template end<0>(); ++el_it) {
		const auto & gtype = el_it->geometry().type();
		const auto & rule = Dune::QuadratureRules<double, dim>::rule(gtype,
				order);

		double elemIntegral = 0.0;
		for (auto ii = rule.begin(); ii != rule.end(); ++ii) {
			const auto & xi = ii->position();
			double omegai = ii->weight();
			double detJaci = el_it->geometry().integrationElement(xi);
			Dune::FieldVector<double, 1> functioni = 0.0;
			dgf.evaluate(*el_it, xi, functioni);
			elemIntegral += functioni * omegai * detJaci;
		}
		integral += elemIntegral;
	}
	return integral;
}

template<typename DGF>
double Integration::boundary_integral(DGF const & dgf, int order) {
	double integral = 0.0;
	auto const & gv = dgf.getGridView();
	const int dim = gv.dimension;

	// By all elements
	auto el_it = gv.template begin<0>();
	for (; el_it != gv.template end<0>(); ++el_it) {
		const auto el_geo = el_it->geometry();
		auto isit_end = gv.iend(*el_it);
		auto isit = gv.ibegin(*el_it);

		// By all sides
		for (; isit != isit_end; ++isit) {
			if (isit->boundary()) {
				auto outerNormal = isit->centerUnitOuterNormal();
				const auto igeo = isit->geometry();  // intersection geometry
				const auto gt = igeo.type();
				const auto& rule = Dune::QuadratureRules<double, dim - 1>::rule(
						gt, order);

				double side_integral = 0.0;

				auto iq = rule.begin();
				for (; iq != rule.end(); ++iq) {
					typename DGF::Traits::RangeType fval;
					dgf.evaluate(*el_it,
							el_geo.local(igeo.global(iq->position())), fval);
					double weight = iq->weight();
					// | det (grad g) |
					double detjac = igeo.integrationElement(iq->position());
					side_integral += (fval * outerNormal) * weight * detjac;
				}
				integral += side_integral;
				//std::cout << "integral = " << integral << std::endl;
			}
		}  // end by all sides
	}  // end by all elements
	return integral;
}



#endif /* SRC_INTEGRATION_HH_ */
