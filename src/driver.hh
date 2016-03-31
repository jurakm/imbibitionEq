#ifndef _DRIVER_HH_17294361_
#define _DRIVER_HH_17294361_

#include <chrono>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
//#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/common/constraints.hh>  // constraints function
#include <dune/pdelab/gridfunctionspace/interpolate.hh>   // interpolate function

#include "parameters.hh"
#include "integration.hh"
#include "bctype.hh"
#include "space_operator.hh"
#include "time_operator.hh"
#include "test.hh"
#include "text_output.hh"
#include "timeMng.hh"
#include "linear_analytic.hh"

/** \brief Driver routine coordinating all calculations except the grid construction.
 *
 *  @tparam GV = Leaf grid view tip 
 *  @tparam Params = Parameters class
 *
 *  @param gv = leaf grid view 
 *  @param params = problem parameters
 *  */
template<typename GV, typename Params>
int driver(GV const& gv, Params params)  // take a copy of params
{
	// analytic solution goes to special driver
	if(params.model == Params::analytic_const || params.model == Params::analytic_var
			                                  || params.model == Params::analytic_new
			                                  ){
		return lin_analytic_driver(params);
	}
	// <<<1>>>
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;

	auto start = std::chrono::system_clock::now();
	TimeMng<Real> timeMng(params);

	// <<<2>>> Grid function space
	constexpr int fem_order = 1;  // 1 or 2 otherwise not implemented
	typedef Dune::PDELab::QkLocalFiniteElementMap<GV, Coord, Real, fem_order> FEM;
	FEM fem(gv);
	typedef Dune::PDELab::ConformingDirichletConstraints CON;
	typedef Dune::PDELab::ISTLVectorBackend<> VBE;
//	typedef Dune::PDELab::istl::VectorBackend<> VBE;
	typedef Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE> GFS;
	GFS gfs(gv, fem);

	BCTypeParam bctype;
	bctype.setTime(timeMng.time);
	typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
	CC cc;
	Dune::PDELab::constraints(bctype, gfs, cc);

	// <<<3>>> Instationary space local operator

	typedef SpaceLocalOperator<BCTypeParam, Params> SLOP;
	SLOP slop(bctype, params, 4);
	// <<<4>>> Time local operator
	typedef TimeLocalOperator<Params> TLOP;
	TLOP tlop(params, 4);
	typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
	MBE mbe(9);
	// <<<5>>>
	typedef Dune::PDELab::GridOperator<GFS, GFS, SLOP, MBE, Real, Real, Real,
			CC, CC> GO0;
	GO0 go0(gfs, cc, gfs, cc, slop, mbe);
	typedef Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, Real, Real, Real,
			CC, CC> GO1;
	GO1 go1(gfs, cc, gfs, cc, tlop, mbe);
	// <<<6>>> I
	typedef Dune::PDELab::OneStepGridOperator<GO0, GO1> IGO;
	IGO igo(go0, go1);

	// <<<7>>> Interpolacija Dirichletovog rubnog te početnog uvjeta
	typedef typename IGO::Traits::Domain U;
//	std::cout<< abi::__cxa_demangle(	typeid(U).name(), 0,0,0) << std::endl;
	U uold(gfs, 0.0);                                       // solution in t=t^n
	BCExtension<GV, Real, Params> g(gv, params);
	g.setTime(timeMng.time);
	Dune::PDELab::interpolate(g, gfs, uold); // the initial and boundary conditions

	// <<<8>>> Select a linear solver backend
	typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
	bool verbose = false; //true;
	LS ls(5000, verbose);

//  // <<<9>>> Solver for linear problem per stage
//  typedef typename IfNonlinear<IGO,LS,U,nonlinear>::PDESOLVER PDESOLVER;
//  typedef Dune::PDELab::StationaryLinearProblemSolver<IGO,LS,U> PDESOLVER;
//  PDESOLVER pdesolver(igo,ls,1e-10);
	// <<<9>>> Solver for non−linear problem per stage
	typedef Dune::PDELab::Newton<IGO, LS, U> PDESOLVER; // In linear and nonlinear case
	PDESOLVER pdesolver(igo, ls);
	pdesolver.setParameters(params.input_data.sub("NewtonParameters"));
//	pdesolver.setReassembleThreshold (0.0);
//	pdesolver.setVerbosityLevel (2);
//	pdesolver.setReduction (1e-10);
//	pdesolver.setMinLinearReduction (1e-4);
//	pdesolver.setMaxIterations (25);
//	pdesolver.setLineSearchMaxIterations (10);

	// <<<10>>> Ovdje se bira vremenska diskretizacija (specijalna RK metoda 2. reda)
//  Dune::PDELab::Alexander2Parameter<Real> method;               // second order
	Dune::PDELab::OneStepThetaParameter<Real> method(1.0); // implicit first order
	Dune::PDELab::OneStepMethod<Real, IGO, PDESOLVER, U, U> osm(method, igo, pdesolver);
	osm.setVerbosityLevel(0);

	std::vector<std::pair<double, double> > volume_values; // (t, int (t))
	std::vector<std::pair<double, double> > bdry_values;  // (t, int (t))

	typedef Dune::PDELab::DiscreteGridFunction<GFS, U> DGF;
	typedef Dune::PDELab::DiscreteGridFunctionGradient<GFS, U> DGFG;

	// File name and variable name depending on the model
	std::string filenm = params.str_sname + params.simulation_names[params.model];
//	if (params.model == Params::nonlinear)
//		filenm = params.str_sname + "nlin";
//	else if (params.model == Params::new_nonlinear)
//		filenm = params.str_sname + "n_nlin";
//	else if (params.model == Params::constant_linear)
//		filenm = params.str_sname + "clin";
//	else if (params.model == Params::variable_linear)
//		filenm = params.str_sname + "vlin";
//	else
//		throw std::runtime_error("Wrong model.");

	Integration integrator;
	// <<<11>>> write out the initial condition
	Dune::PDELab::FilenameHelper fn(filenm);
	{
		DGF udgf(gfs, uold);
		if (params.vtkout) {
			Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, fem_order - 1);
			vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(udgf, filenm));
			vtkwriter.write(fn.getName(), Dune::VTK::ascii);
			fn.increment();
		}
		if (params.txtout) {
			/////////// Gnuplot output
			TextOutput<DGF, GV> textwriter(udgf, gv);
			textwriter.write(filenm + std::to_string(0) + ".txt");
		}
		// Calculate integrals corresponding to the initial condition
		DGFG grad_udgf(gfs, uold);
		integrator.integrate(0.0, udgf, grad_udgf);
//		double intOfsol = volume_integral(udgf, 3);
//		volume_values.push_back(std::make_pair(0.0, intOfsol));
//
//		double intOverBoundary = boundary_integral(grad_udgf, 3);
//		bdry_values.push_back(std::make_pair(0.0, intOverBoundary));
	}

	// <<<12>>> the time loop
	U unew(gfs, 0.0);
	while (!timeMng.done()) {
		++timeMng;
		// postavi novo vrijeme u BC klasu
		bctype.setTime(timeMng.time);
		cc.clear();
		Dune::PDELab::constraints(bctype, gfs, cc);
 //       std::cout << "dt = " << timeMng.dt << " ";
		double dt_applied = osm.apply(timeMng.time - timeMng.dt, timeMng.dt, uold, g, unew);

		if (dt_applied != timeMng.dt) {
			// this probably never happen
			std::cout << "dt changed : old " << timeMng.dt << ", new: "
					<< dt_applied << "\n";
			const double diff = dt_applied - timeMng.dt;
			timeMng.dt = dt_applied;
			timeMng.time += diff;
		}

		//int noIter = osm.getPDESolver().ls_result().iterations;
		int noIter = osm.getPDESolver().result().iterations;
 //       std::cout << "no of nonlinear iterations = " << noIter << "\n";
		// kontrola vremenskog koraka
		timeMng.set_requested_dt(noIter);
		// graphics
		DGF udgf(gfs, unew);
		if (params.vtkout and timeMng.doOutput()) {
			Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, fem_order - 1);
			vtkwriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGF>>(udgf, filenm));
			vtkwriter.write(fn.getName(), Dune::VTK::ascii);
			fn.increment();
		}
		if (params.txtout and timeMng.doOutput()) {
			/////////// Gnuplot output
			const int count = timeMng.output_count;
			TextOutput<DGF, GV> textwriter(udgf, gv);
			textwriter.write(filenm + std::to_string(count) + ".txt");
		}
		////////////////////////////
		if(params.model != Params::new_nonlinear){
		   DGFG grad_udgf(gfs, unew);
		   integrator.integrate(timeMng.time, udgf, grad_udgf);
		}
		else{
			U beta_u(gfs, 0.0);
			auto  bit = beta_u.begin();
			for(auto it = unew.begin(); it != unew.end(); ++it, ++bit) *bit = params.beta(*it);
		    DGFG grad_beta_udgf(gfs, beta_u);
		    integrator.integrate(timeMng.time, udgf, grad_beta_udgf);
 		}

		uold = unew;
		//     std::cout << "t = " << timeMng.time << " (dt = " << timeMng.dt << ")\n";
	}

	integrator.volume_derivative();
	integrator.print(filenm + "-flux.txt", params);
	auto end = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds> (end - start);
	return duration.count();
//	return; // timeMng.output_count;
}

#endif
