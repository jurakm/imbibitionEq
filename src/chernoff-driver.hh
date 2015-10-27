#ifndef _CHERNOFF_DRIVER_HH__
#define _CHERNOFF_DRIVER_HH__

#include <chrono>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
//#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/constraints/common/constraints.hh>  // constraints function
#include <dune/pdelab/gridfunctionspace/interpolate.hh>   // interpolate function

#include "integration.hh"
#include "bctype.hh"
#include "chernoff_operator.hh"
#include "test.hh"
#include "text_output.hh"
#include "timeMng.hh"
#include "linear_analytic.hh"

/** Driver routine coordinating all calculations except the grid construction.
 *
 *  @tparam GV = Leaf grid view tip 
 *  @tparam Params = Parameters class
 *
 *  @param gv = leaf grid view 
 *  @param params = problem parameters
 *  */
template<typename GV, typename Params>
int chernoff_driver(GV const& gv, Params params)  // take a copy of params
{
    // analytic solution goes to special driver
    if(params.model == Params::analytic_const || params.model == Params::analytic_var
                                              || params.model == Params::analytic_new
                                              || params.model == Params::analytic_new1){
        return lin_analytic_driver(params);
    }
    // <<<1>>>
    typedef typename GV::Grid::ctype Coord;
    typedef double Real;
    constexpr int fem_order = 1;  // 1 or 2 otherwise not implemented
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV, Coord, Real, fem_order> FEM;
    typedef Dune::PDELab::ConformingDirichletConstraints CON;
    typedef Dune::PDELab::ISTLVectorBackend<> VBE;
    typedef Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE> GFS;
    typedef typename Dune::PDELab::BackendVectorSelector<GFS,Real>::Type U;
    typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
    typedef Dune::PDELab::DiscreteGridFunction<GFS, U> DGF;
    typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS, U> DGFG;
    typedef LocalOperator<BCTypeParam, Params, DGF> LOP;
    typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
    typedef Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, Real, Real, Real, CC, CC> GO;
    typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
    typedef Dune::PDELab::Newton<GO, LS, U> PDESOLVER; // In linear and nonlinear case
    
    auto start = std::chrono::system_clock::now();
    TimeMng<Real> timeMng(params);

    FEM fem(gv);
    GFS gfs(gv, fem);
    BCTypeParam bctype;
    bctype.setTime(timeMng.time);
    CC cc;
    Dune::PDELab::constraints(bctype, gfs, cc);

//  std::cout<< abi::__cxa_demangle(    typeid(U).name(), 0,0,0) << std::endl;
    U uold(gfs, 0.0);                                       // solution in t=t^n
    BCExtension<GV, Real, Params> g(gv, params);
    g.setTime(timeMng.time);
    Dune::PDELab::interpolate(g, gfs, uold); // the initial and boundary conditions
    DGF udgf_old(gfs, uold);
    
    U unew(gfs, 0.0);
    DGF udgf_new(gfs, unew);
    DGFG grad_udgf_new(gfs, unew);
    
    MBE mbe(9);
    LOP lop(bctype, params, udgf_old, timeMng.dt, 4);
    GO go(gfs, cc, gfs, cc, lop, mbe);

    bool verbose = false; //true;
    LS ls(5000, verbose);
    PDESOLVER pdesolver(go, unew, ls);
    pdesolver.setParameters(params.input_data.sub("NewtonParameters"));
    pdesolver.setVerbosityLevel(2);

    std::vector<std::pair<double, double> > volume_values; // (t, int (t))
    std::vector<std::pair<double, double> > bdry_values;  // (t, int (t))
 
 //   typename GO::Traits::Jacobian jac(go);
 //   std::cout << jac.patternStatistics()<< std::endl;


    // File name and variable name depending on the model
    std::string filenm = params.str_sname + params.simulation_names[params.model];
//  if (params.model == Params::nonlinear)
//      filenm = params.str_sname + "nlin";
//  else if (params.model == Params::new_nonlinear)
//      filenm = params.str_sname + "n_nlin";
//  else if (params.model == Params::constant_linear)
//      filenm = params.str_sname + "clin";
//  else if (params.model == Params::variable_linear)
//      filenm = params.str_sname + "vlin";
//  else
//      throw std::runtime_error("Wrong model.");

    Integration integrator;
    // <<<11>>> write out the initial condition
    Dune::PDELab::FilenameHelper fn(filenm);
    {
        if (params.vtkout) {
            Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, fem_order - 1);
            vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(udgf_old, filenm));
            vtkwriter.write(fn.getName(), Dune::VTK::ascii);
            fn.increment();
        }
        if (params.txtout) {
            /////////// Gnuplot output
            TextOutput<DGF, GV> textwriter(udgf_old, gv);
            textwriter.write(filenm + std::to_string(0) + ".txt");
        }
        // Calculate integrals corresponding to the initial condition
        DGFG grad_udgf_old(gfs, uold);
        integrator.integrate(0.0, udgf_old, grad_udgf_old);
//      double intOfsol = volume_integral(udgf, 3);
//      volume_values.push_back(std::make_pair(0.0, intOfsol));
//
//      double intOverBoundary = boundary_integral(grad_udgf, 3);
//      bdry_values.push_back(std::make_pair(0.0, intOverBoundary));
    }

    // <<<12>>> the time loop
  
    
    while (!timeMng.done()) {
        ++timeMng;
        // postavi novo vrijeme u BC klasu
        bctype.setTime(timeMng.time);
        lop.setDt(timeMng.dt);
        lop.setTime(timeMng.time);
//        cc.clear();
//        Dune::PDELab::constraints(bctype, gfs, cc);
        pdesolver.apply();

        int noIter = pdesolver.result().iterations;
        std::cout << "no of nonlinear iterations = " << noIter << "\n";
        // kontrola vremenskog koraka
        timeMng.set_requested_dt(noIter);
        // graphics
        if (params.vtkout and timeMng.doOutput()) {
            Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, fem_order - 1);
            vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(udgf_new, filenm));
            vtkwriter.write(fn.getName(), Dune::VTK::ascii);
            fn.increment();
        }
        if (params.txtout and timeMng.doOutput()) {
            /////////// Gnuplot output
            const int count = timeMng.output_count;
            TextOutput<DGF, GV> textwriter(udgf_new, gv);
            textwriter.write(filenm + std::to_string(count) + ".txt");
        }
        ////////////////////////////
        if(params.model != Params::new_nonlinear){
           integrator.integrate(timeMng.time, udgf_new, grad_udgf_new);
        }
        else{
            U beta_u(gfs, 0.0);
            auto  bit = beta_u.begin();
            for(auto it = unew.begin(); it != unew.end(); ++it, ++bit) *bit = params.beta(*it);
            DGFG grad_beta_udgf(gfs, beta_u);
            integrator.integrate(timeMng.time, udgf_new, grad_beta_udgf);
        }

        uold = unew;
        //     std::cout << "t = " << timeMng.time << " (dt = " << timeMng.dt << ")\n";
    }

    integrator.volume_derivative();
    integrator.print(filenm + "-flux.txt", params);
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds> (end - start);
    return duration.count();
//  return; // timeMng.output_count;
}

#endif
