#ifndef _CHERNOFF_OPERATOR__
#define _CHERNOFF_OPERATOR__

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/idefault.hh>

#include <dune/pdelab/common/geometrywrapper.hh>
//#include <dune/grid/common/geometry.hh>

#include <stdexcept>
#include <cassert>

#include "parameters.hh"
/** @brief Local operator for Chernoff formula method
 * 
 * A local operator for solving the equation 
 *
 *  \f[ \Phi \partial S/\partial t - \delta^2 k \Delta \beta(S) = 0  \quad \mbox{in }\; Y \times  (0,T)\f]
 *    \f[               S(x,t) = g(t)   \quad \mbox{on }\; \partial Y\times (0,T) \f]
 *    \f[               S(x,0) = g(0)   \quad \mbox{on }\; Y\f]
 * 
 * with conforming finite elements. We use backward Euler discrtization leading to the equation:
 * 
 * \f[ \Phi \frac{S^n - S^{n-1}}{\Delta t}  - \delta^2 k \Delta \beta(S^n) = 0 \quad \mbox{in }\; Y \times  (0,T)\f]
 *    \f[               S^n = g(t^n)   \quad \mbox{on }\; \partial Y \times (0,T) \f]
 *  for \f$n=1,2,\ldots\f$, where \f$S^0=g(0)\f$. 
 * 
 *
 * \tparam BCType parameter class indicating the type of boundary condition
 */
template<typename BCType, typename Params, typename DGF>
class LocalOperator: 
        public Dune::PDELab::NumericalJacobianApplyVolume< LocalOperator<BCType, Params, DGF> >,
        public Dune::PDELab::NumericalJacobianVolume<LocalOperator<BCType, Params, DGF> >,
        public Dune::PDELab::NumericalJacobianApplyBoundary<LocalOperator<BCType, Params, DGF> >,
        public Dune::PDELab::NumericalJacobianBoundary<LocalOperator<BCType, Params, DGF> >,
        public Dune::PDELab::FullVolumePattern,
        public Dune::PDELab::LocalOperatorDefaultFlags {
public:
    enum { doPatternVolume = true };
    enum { doAlphaVolume = true };
    enum { doLambdaVolume = true };

    /// Constructor
    LocalOperator(const BCType& bctype_, const Params& coeff_, DGF const & dgf, double  dt, unsigned int intorder_ = 3) :
            time_(0.0), dt_(dt), bctype(bctype_), coeff(coeff_), dgf_(dgf),  intorder(intorder_) {
    }
    
    /** Set dt for new calculation. */
    void setDt(double dt) { dt_ = dt; assert(dt_ > 0.0); }
    
    /** Set time for subsequent calculation. */
    void setTime(double time) { time_ = time; }
    
    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const {
        // assume Galerkin: lfsu == lfsv
        // This yields more efficient code since the local function space only
        // needs to be evaluated once, but would be incorrect for a finite volume
        // method

        // dimensions
        const int dim  = EG::Geometry::mydimension;
        const int dimw = EG::Geometry::coorddimension;

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

        const double poro = coeff.poro;
        // loop over quadrature points
        for (auto it = rule.begin(); it != rule.end(); ++it) {
            // evaluate basis functions on reference element
            std::vector<Range> phi(lfsu.size());
            lfsu.finiteElement().localBasis().evaluateFunction(it->position(), phi);

            // compute solution u (or beta(u)) at integration point
            RF u = 0.0;
            for (size_type i = 0; i < lfsu.size(); ++i)
                    u += x(lfsu, i) * phi[i];

            double alpha = coeff.k * coeff.delta * coeff.delta;

            if (coeff.model == Params::nonlinear)
                alpha *= coeff.alpha(u);  // ovo nije moguće
            else
                alpha *= coeff.a_g(time_);

//          else if (coeff.model == Params::constant_linear)
//              alpha *= coeff.mean_alpha;
//          else if (coeff.model == Params::variable_linear)
//              alpha *= coeff.alpha(coeff.bdry(time_));
//          else {
//              // if new_nonlinear alpha is not changed
//          }
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
            if (coeff.model == Params::new_nonlinear || coeff.model == Params::Params::chernoff) {
                for (size_type i = 0; i < lfsu.size(); ++i)
                    gradu.axpy(coeff.beta(x(lfsu, i)), gradphi[i]); // grad beta(u)
            } else {  // In all other cases calculate grad u.
                for (size_type i = 0; i < lfsu.size(); ++i)
                    gradu.axpy(x(lfsu, i), gradphi[i]);
            }
            // integrate grad u * grad phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i = 0; i < lfsu.size(); ++i)
                r.accumulate(lfsu, i, poro * ( u/dt_) * phi[i] * factor  +  alpha * (gradu * gradphi[i]) * factor);
        }
    }
     // volume integral depending only on test functions
     // Access to dgf is localized here.
     template<typename EG, typename LFSV, typename R>
     void lambda_volume ( const EG& eg, const LFSV& lfsv, R& r ) const {
       // dimensions
        const int dim  = EG::Geometry::mydimension;
        const int dimw = EG::Geometry::coorddimension;

        // extract some types
        typedef typename LFSV::Traits::FiniteElementType::Traits::LocalBasisType::Traits LBTraits;
        typedef typename LBTraits::DomainFieldType DF;
        typedef typename LBTraits::RangeFieldType  RF;
        typedef typename LBTraits::RangeType       Range;
        typedef typename LFSV::Traits::SizeType   size_type;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF, dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt, intorder);

        const double poro = coeff.poro;
        // loop over quadrature points
        for (auto it = rule.begin(); it != rule.end(); ++it) {
            // evaluate basis functions on reference element
            std::vector<Range> phi(lfsv.size());
            lfsv.finiteElement().localBasis().evaluateFunction(it->position(), phi);

            typename DGF::Traits::RangeType u_old;
            dgf_.evaluate(eg.entity(), it->position(), u_old);
            // integrate grad u * grad phi_i
            RF factor = it->weight() * eg.geometry().integrationElement(it->position());
            for (size_type i = 0; i < lfsv.size(); ++i)
                r.accumulate(lfsv, i, - poro * ( u_old/dt_) * phi[i] * factor );
        }
     }

protected:
    double time_;
    double dt_;

private:
    const BCType& bctype;
    const Params & coeff;
    /** Discrete grid function representing the solution at preceeding time level. */
    DGF const & dgf_;     
    unsigned int intorder;
};



#endif