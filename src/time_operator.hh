#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/idefault.hh>
//#include <dune/grid/common/geometry.hh>

/** \brief Bilinear form under the time derivative.
 *
 *
 * Class representing local operator for the accumulation term
 * \f[
         \int_\Omega \Phi u v dx
 * \f]
 * where \f$\Phi\f$ is the porosity.
 *
 * \tparam Params Parameter class here used only to get the porosity.
 */
template <typename Params>
class TimeLocalOperator 
  : public Dune::PDELab::NumericalJacobianApplyVolume<TimeLocalOperator<Params> >,
    public Dune::PDELab::NumericalJacobianVolume<TimeLocalOperator<Params> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
public:
  //! pattern assembly flags
  enum { doPatternVolume = true };

  //! residual assembly flags
  enum { doAlphaVolume = true };

  /// Constructor.
  TimeLocalOperator (Params const & params, unsigned int intorder_=2)
    : params_(params), intorder(intorder_), time(0.0)
  {}

  //! set time for subsequent evaluation
  void setTime (double t) {time = t;}

  //! volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // domain and range field type
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
    typedef typename LFSU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
    typedef typename LFSU::Traits::SizeType size_type;
        
    // dimensions
    const int dim = EG::Geometry::mydimension;
    const double poro = params_.poro;
    // select quadrature rule
    Dune::GeometryType gt = eg.geometry().type();
    const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,intorder);

    // loop over quadrature points
    for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
      {
        // evaluate basis functions
        std::vector<RangeType> phi(lfsu.size());
        lfsu.finiteElement().localBasis().evaluateFunction(it->position(),phi);

        // evaluate u
        RF u=0.0;
        for (size_type i=0; i<lfsu.size(); i++)
          u += x(lfsu,i)*phi[i];

        // u*phi_i
        RF dx = it->weight() * eg.geometry().integrationElement(it->position());
        for (size_type i=0; i<lfsu.size(); i++)
          r.accumulate(lfsu, i, poro * u * phi[i] * dx);
      }
  }
private:
  Params params_;
  unsigned int intorder;
  double time;
};
