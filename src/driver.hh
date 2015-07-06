#ifndef _DRIVER_HH_17294361_
#define _DRIVER_HH_17294361_

#include <dune/pdelab/finiteelementmap/qkfem.hh>
//#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/newton/newton.hh>

#include "integration.hh"
#include "bctype.hh"
#include "space_operator.hh"
#include "time_operator.hh"
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
void driver(GV const& gv, Params params)  // take a copy of params
{
	// analytic solution goes to special driver
	if(params.model == Params::analytic_const){
		lin_analytic(params);
		return;
	}
	// <<<1>>>
	typedef typename GV::Grid::ctype Coord;
	typedef double Real;

	TimeMng<Real> timeMng(params);

	// <<<2>>> Grid function space
	constexpr int fem_order = 1;  // 1 or 2 otherwise not implemented
	typedef Dune::PDELab::QkLocalFiniteElementMap<GV, Coord, Real, fem_order> FEM;
	FEM fem(gv);
	typedef Dune::PDELab::ConformingDirichletConstraints CON;
	typedef Dune::PDELab::ISTLVectorBackend<> VBE;
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
			vtkwriter.addVertexData(
					new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf, filenm));
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
//      std::cout << "no of nonlinear iterations = " << noIter << "\n";
		// kontrola vremenskog koraka
		timeMng.set_requested_dt(noIter);
		// graphics
		DGF udgf(gfs, unew);
		DGFG grad_udgf(gfs, unew);
		if (params.vtkout and timeMng.doOutput()) {
			Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, fem_order - 1);
			vtkwriter.addVertexData(
					new Dune::PDELab::VTKGridFunctionAdapter<DGF>(udgf, filenm));
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
		integrator.integrate(timeMng.time, udgf, grad_udgf);

//		double intOfsol = volume_integral(udgf, 3);
//		volume_values.push_back(std::make_pair(timeMng.time, intOfsol));
////      std::cout << "(t, int(u)) = (" << timeMng.time <<"," << std::setprecision(12) << intOfsol << ")\n";
//		double intOverBoundary = boundary_integral(grad_udgf, 3);
//		bdry_values.push_back(std::make_pair(timeMng.time, intOverBoundary));
////      std::cout << "(t, int_bdry(grad u)) = (" << timeMng.time <<"," << std::setprecision(12) << intOverBoundary << ")\n";

		uold = unew;
		//     std::cout << "t = " << timeMng.time << " (dt = " << timeMng.dt << ")\n";
	}

	integrator.volume_derivative();
//	std::vector<std::pair<double, double> > volume_values_der(
//			volume_values.size()); // (t, d/dt int S(t))
//
//	double dt0 = volume_values[1].first - volume_values[0].first;
//	double dS0 = volume_values[1].second - volume_values[0].second;
//	volume_values_der[0] = std::make_pair(volume_values[0].first, dS0 / dt0);
//
//	unsigned int nn = volume_values.size() - 1;
//	for (unsigned int i = 1; i < nn; ++i) {
//		double dt = volume_values[i + 1].first - volume_values[i - 1].first;
//		double dS = volume_values[i + 1].second - volume_values[i - 1].second;
//
//		volume_values_der[i] = std::make_pair(volume_values[i].first, dS / dt);
//	}
//	double dtn = volume_values[nn].first - volume_values[nn - 1].first;
//	double dSn = volume_values[nn].second - volume_values[nn - 1].second;
//	volume_values_der[nn] = std::make_pair(volume_values[nn].first, dSn / dtn);


	integrator.print(filenm + "-flux.txt", params);
//	std::ofstream out1(filenm + "-flux.txt");
//
//	out1
//			<< "#     t        Phi d/dt int S   k delta^2 alpha(g(t)) int bdry grad S .n    bdry(t)\n";
//	const double kdd = params.k * params.delta * params.delta;
//	for (unsigned int i = 0; i < volume_values.size(); ++i) {
//		double time = volume_values_der[i].first;
//		if (std::abs(time - bdry_values[i].first) >= 1.0E-12)
//			throw std::runtime_error(
//					std::string("Time error! i = ") + std::to_string(i) + " "
//							+ std::to_string(time) + " "
//							+ std::to_string(bdry_values[i].first));
//
//		out1 << std::setw(10) << std::setprecision(4) << time
//				<< " "  // time
//				<< std::setw(12) << std::setprecision(6)
//				<< params.poro * volume_values_der[i].second; //  " = d/dt int S"
//		if (params.model == Params::nonlinear)
//			out1 << "                 " << std::setw(12) << std::setprecision(6)
//					<< params.alpha(params.bdry(time)) * kdd
//							* bdry_values[i].second;
//		// "= k delta^2 alpha(g(t)) int bdry grad S .n  "
//		else if (params.model == Params::variable_linear)
//			out1 << "                 " << std::setw(12) << std::setprecision(6)
//					<< params.alpha_reg(params.bdry(time)) * kdd
//							* bdry_values[i].second;
//		// "= k delta^2 alpha(g(t)) int bdry grad S .n  "
//		else if (params.model == Params::constant_linear)
//			out1 << "                 " << std::setw(12) << std::setprecision(6)
//					<< params.mean_alpha * kdd * bdry_values[i].second;
//		// "= k delta^2 mean_alpha int bdry grad S .n  "
//
//		out1 << "             " << std::setw(12) << std::setprecision(6)
//				<< params.bdry(time) << "\n";
//	}
//	out1.close();

	return; // timeMng.output_count;
}

#endif
