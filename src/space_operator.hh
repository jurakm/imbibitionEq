#include <dune/pdelab/localoperator/idefault.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>

#include <stdexcept>
#include "parameters.hh"
/** Space part of a local operator for solving the equation
 *
 *   Phi dS/dt - delta^2 k div(a(S) grad S) = 0   in \Omega x (0,T)
 *                  S = g   on \partial\Omega x (0,T)
 *
 * with conforming finite elements
 *
 * \tparam BCType parameter class indicating the type of boundary condition
 */
template<class BCType, class Params>
class StationaryLocalOperator: public Dune::PDELab::NumericalJacobianApplyVolume<
		StationaryLocalOperator<BCType, Params> >,
		public Dune::PDELab::NumericalJacobianVolume<StationaryLocalOperator<BCType, Params> >,
		public Dune::PDELab::NumericalJacobianApplyBoundary<StationaryLocalOperator<BCType, Params> >,
		public Dune::PDELab::NumericalJacobianBoundary<StationaryLocalOperator<BCType, Params> >,
		public Dune::PDELab::FullVolumePattern,
		public Dune::PDELab::LocalOperatorDefaultFlags {
public:
	// pattern assembly flags
	enum {
		doPatternVolume = true
	};

	// residual assembly flags
	enum {
		doAlphaVolume = true
	};
    /// Constructor
	StationaryLocalOperator(const BCType& bctype_,
			const Params& coeff_,
			unsigned int intorder_ = 3) :
			time_(0.0), bctype(bctype_), coeff(coeff_), intorder(intorder_) {
	}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x,
			const LFSV& lfsv, R& r) const {
		// assume Galerkin: lfsu == lfsv
		// This yields more efficient code since the local function space only
		// needs to be evaluated once, but would be incorrect for a finite volume
		// method

		// dimensions
		const int dim  = EG::Geometry::dimension;
		const int dimw = EG::Geometry::dimensionworld;

		// extract some types
		typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;
		typedef typename LBTraits::DomainFieldType DF;
		typedef typename LBTraits::RangeFieldType  RF;
		typedef typename LBTraits::JacobianType    Jacobian;
		typedef typename LBTraits::RangeType       Range;
		typedef Dune::FieldVector<RF, dimw>        Gradient;
		typedef typename LFSU::Traits::SizeType    size_type;

		// select quadrature rule
		Dune::GeometryType gt = eg.geometry().type();
		const Dune::QuadratureRule<DF, dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt, intorder);

		// loop over quadrature points
		for (auto it = rule.begin(); it != rule.end(); ++it) {
			// evaluate basis functions on reference element
			std::vector<Range> phi(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateFunction(it->position(), phi);

			// compute solution u (or beta(u)) at integration point
			RF u = 0.0;
			if (coeff.model == Params::nonlinear) {
				for (size_type i = 0; i < lfsu.size(); ++i)
					u += x(lfsu, i) * phi[i];
			}

			double alpha = coeff.k * coeff.delta * coeff.delta;

			if (coeff.model == Params::nonlinear)
				alpha *= coeff.alpha(u);
			else if (coeff.model == Params::constant_linear)
				alpha *= coeff.mean_alpha;
			else if (coeff.model == Params::variable_linear)
				alpha *= coeff.alpha(coeff.bdry(time_));
			else {
				// if new_nonlinear alpha is not changed
			}
			// evaluate gradient of basis functions on reference element
			std::vector<Jacobian> js(lfsu.size());
			lfsu.finiteElement().localBasis().evaluateJacobian(it->position(), js);

			// transform gradients from reference element to real element
			const Dune::FieldMatrix<DF, dimw, dim> &jac = eg.geometry().jacobianInverseTransposed(it->position());
			std::vector<Gradient> gradphi(lfsu.size());
			for (size_type i = 0; i < lfsu.size(); i++)
				jac.mv(js[i][0], gradphi[i]);

			// compute gradient of u or gradient of beta(u)
			Gradient gradu(0.0);
			if (coeff.model == Params::new_nonlinear) {
				for (size_type i = 0; i < lfsu.size(); ++i)
					gradu.axpy(coeff.beta(x(lfsu, i)), gradphi[i]); // grad beta(u)
			} else {  // In all other cases calculate grad u.
				for (size_type i = 0; i < lfsu.size(); ++i)
					gradu.axpy(x(lfsu, i), gradphi[i]);
			}
			// integrate grad u * grad phi_i
			RF factor = it->weight() * eg.geometry().integrationElement(it->position());
			for (size_type i = 0; i < lfsu.size(); ++i)
				r.accumulate(lfsu, i, alpha * (gradu * gradphi[i]) * factor);
		}
	}

protected:
	double time_;

private:
	const BCType& bctype;
	const Params & coeff;
	unsigned int intorder;
};

/** Space part of a local operator for solving the equation
 *
 *   Phi dS/dt - delta^2 k div(a(S) grad S) = 0   in \Omega x (0,T)
 *                  S = g   on \partial\Omega x (0,T)
 *
 * with conforming finite elements
 *
 *
 * \tparam B = class indicating the type of boundary condition
 * \tparam C =  class of problem parameters
 */
template<class BCType, class Params>
class SpaceLocalOperator: public StationaryLocalOperator<BCType, Params>,
		public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double> // default methods
{
	BCType& b;
public:
	SpaceLocalOperator(BCType& b_, const Params& c_, unsigned int intorder_ = 2) :
			StationaryLocalOperator<BCType, Params>(b_, c_, intorder_), b(b_) {
	}

	void preStep(double time, double dt, int stages) {
		this->time_ = time;
//    b.setTime(time); // postavi korektno vrijeme rubnom uvjetu
		Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>::preStep( time, dt, stages);
	}
};
