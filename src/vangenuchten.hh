#ifndef _TWOPHASEPARAM_HH_IS_INCLUDED_
#define _TWOPHASEPARAM_HH_IS_INCLUDED_

#include <cassert>
#include <cmath>


namespace Dumux
{
class VanGenuchtenParams
{
public:
    using Scalar = double;
    VanGenuchtenParams(){}
    VanGenuchtenParams(Scalar vgAlpha, Scalar vgn) : vgAlpha_(vgAlpha), vgn_(vgn), vgm_(1 -1.0/vgn)
    {}
    /*!
     * \brief Return the \f$\alpha\f$ shape parameter [1/Pa] of van Genuchten's
     *        curve.
     */
    Scalar vgAlpha() const
    { return vgAlpha_; }

    /*!
     * \brief Set the \f$\alpha\f$ shape parameter [1/Pa] of van Genuchten's
     *        curve.
     */
    void setVgAlpha(Scalar v)
    { vgAlpha_ = v; }

    /*!
     * \brief Return the \f$m\f$ shape parameter [-] of van Genuchten's
     *        curve.
     */
    Scalar vgm() const
    { return vgm_; }

    /*!
     * \brief Set the \f$m\f$ shape parameter [-] of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$n = \frac{1}{1 - m}\f$
     */
    void setVgm(Scalar m)
    { vgm_ = m; vgn_ = 1/(1 - vgm_); }

    /*!
     * \brief Return the \f$n\f$ shape parameter [-] of van Genuchten's
     *        curve.
     */
    Scalar vgn() const
    { return vgn_; }

    /*!
     * \brief Set the \f$n\f$ shape parameter [-] of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$m = 1 - \frac{1}{n}\f$
     */
    void setVgn(Scalar n)
    { vgn_ = n; vgm_ = 1 - 1/vgn_; }

private:
    Scalar vgAlpha_ = 0.0;
    Scalar vgn_ = 0.0;
    Scalar vgm_ = 0.0;
};

using std::pow;

class VanGenuchten{
public:
  using Params =  VanGenuchtenParams;
  using Scalar = double;
  /*!
      * \brief The capillary pressure-saturation curve according to van Genuchten.
      *
      * Van Genuchten's empirical capillary pressure <-> saturation
      * function is given by
      * \f[
      p_C = (\overline{S}_w^{-1/m} - 1)^{1/n}/\alpha
      \f]
      * \param swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$
      * \param params    A container object that is populated with the appropriate coefficients for the respective law.
      *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
      *                  is constructed accordingly. Afterwards the values are set there, too.
      */
     static Scalar pc(const Params &params, Scalar swe)
     {
         assert(0 <= swe && swe <= 1);
         return pow(pow(swe, -1.0/params.vgm()) - 1, 1.0/params.vgn())/params.vgAlpha();
     }

     /*!
      * \brief The saturation-capillary pressure curve according to van Genuchten.
      *
      * This is the inverse of the capillary pressure-saturation curve:
      * \f[
      \overline{S}_w = {p_C}^{-1} = ((\alpha p_C)^n + 1)^{-m}
      \f]
      *
      * \param pc        Capillary pressure \f$p_C\f$
      * \param params    A container object that is populated with the appropriate coefficients for the respective law.
      *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
      *                  is constructed accordingly. Afterwards the values are set there, too.
      * \return          The effective saturation of the wetting phase \f$\overline{S}_w\f$
      */
     static Scalar sw(const Params &params, Scalar pc)
     {
         assert(pc >= 0);

         return pow(pow(params.vgAlpha()*pc, params.vgn()) + 1, -params.vgm());
     }

     /*!
      * \brief The partial derivative of the capillary
      *        pressure w.r.t. the effective saturation according to van Genuchten.
      *
      * This is equivalent to
      * \f[
      \frac{\partial p_C}{\partial \overline{S}_w} =
      -\frac{1}{\alpha} (\overline{S}_w^{-1/m} - 1)^{1/n - }
      \overline{S}_w^{-1/m} / \overline{S}_w / m
      \f]
      *
      * \param swe       Effective saturation of the wetting phase \f$\overline{S}_w\f$
      * \param params    A container object that is populated with the appropriate coefficients for the respective law.
      *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
      *                  is constructed accordingly. Afterwards the values are set there, too.
     */
     static Scalar dpc_dsw(const Params &params, Scalar swe)
     {
         assert(0 <= swe && swe <= 1);

         Scalar powSwe = pow(swe, -1/params.vgm());
         return - 1.0/params.vgAlpha() * pow(powSwe - 1, 1.0/params.vgn() - 1)/params.vgn()
             * powSwe/swe/params.vgm();
     }

     /*!
      * \brief The partial derivative of the effective
      *        saturation to the capillary pressure according to van Genuchten.
      *
      * \param pc        Capillary pressure \f$p_C\f$
      * \param params    A container object that is populated with the appropriate coefficients for the respective law.
      *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
      *                  is constructed accordingly. Afterwards the values are set there, too.
      */
     static Scalar dsw_dpc(const Params &params, Scalar pc)
     {
         assert(pc >= 0);

         Scalar powAlphaPc = pow(params.vgAlpha()*pc, params.vgn());
         return -pow(powAlphaPc + 1, -params.vgm()-1)*
             params.vgm()*powAlphaPc/pc*params.vgn();
     }

     /*!
      * \brief The relative permeability for the wetting phase of
      *        the medium implied by van Genuchten's
      *        parameterization.
      *
      * \param swe        The mobile saturation of the wetting phase.
      * \param params    A container object that is populated with the appropriate coefficients for the respective law.
      *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
      *                  is constructed accordingly. Afterwards the values are set there, too.     */
     static Scalar krw(const Params &params, Scalar swe)
     {
         assert(0 <= swe && swe <= 1);

         Scalar r = 1.0 - pow(1.0 - pow(swe, 1.0/params.vgm()), params.vgm());
         return sqrt(swe)*r*r;
     }

     /*!
      * \brief The derivative of the relative permeability for the
      *        wetting phase in regard to the wetting saturation of the
      *        medium implied by the van Genuchten parameterization.
      *
      * \param swe       The mobile saturation of the wetting phase.
      * \param params    A container object that is populated with the appropriate coefficients for the respective law.
      *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
      *                  is constructed accordingly. Afterwards the values are set there, too.
      */
     static Scalar dkrw_dsw(const Params &params, Scalar swe)
     {
         assert(0 <= swe && swe <= 1);

         const Scalar x = 1.0 - std::pow(swe, 1.0/params.vgm());
         const Scalar xToM = std::pow(x, params.vgm());
         return (1.0 - xToM)/std::sqrt(swe) * ( (1.0 - xToM)/2 + 2*xToM*(1.0-x)/x );
     }

     /*!
      * \brief The relative permeability for the non-wetting phase
      *        of the medium implied by van Genuchten's
      *        parameterization.
      *
      * \param swe        The mobile saturation of the wetting phase.
      * \param params    A container object that is populated with the appropriate coefficients for the respective law.
      *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
      *                  is constructed accordingly. Afterwards the values are set there, too.
      */
     static Scalar krn(const Params &params, Scalar swe)
     {
         assert(0 <= swe && swe <= 1);

         return
             pow(1 - swe, 1.0/3) *
             pow(1 - pow(swe, 1.0/params.vgm()), 2*params.vgm());
     }

     /*!
      * \brief The derivative of the relative permeability for the
      *        non-wetting phase in regard to the wetting saturation of
      *        the medium as implied by the van Genuchten
      *        parameterization.
      *
      * \param swe        The mobile saturation of the wetting phase.
      * \param params    A container object that is populated with the appropriate coefficients for the respective law.
      *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
      *                  is constructed accordingly. Afterwards the values are set there, too.
      */
     static Scalar dkrn_dsw(const Params &params, Scalar swe)
     {
         assert(0 <= swe && swe <= 1);

         const Scalar x = std::pow(swe, 1.0/params.vgm());
         return
             -std::pow(1.0 - x, 2*params.vgm())
             *std::pow(1.0 - swe, -2.0/3)
             *(1.0/3 + 2*x/swe);
     }

};
}
#endif
