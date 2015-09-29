// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

 \mainpage Solution of nonlinear imbibition equation and two of its linearizations.

 \section Problem

   We solve so called imbibition equation and we calculate first the solution of the nonlinear
   problem:
   \f[
       \Phi \frac{\partial S}{\partial t} - \delta^2 k_m \Delta \beta(S) = 0 \quad {\rm in }\;  Y,\; t>0
   \f]
   \f[  S = {\cal P}(S^f)\quad{\rm on } \;\partial Y,\; t > 0 \f]
   \f[  S = S_0 \quad {\rm on}\; Y,\; t = 0 \f]

 */
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#define O2SCL_USE_GSL_HANDLER

#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
//#include <dune/common/static_assert.hh>
#include <dune/common/timer.hh>
// For reading of an .input file.
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
//#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
//#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istl.hh>
//#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include "parameters.hh"
//#include "read_input.hh"
#include "driver.hh"
#include "mgf.hh"
#include "linear_analytic.hh"
#include "imbibitionFunctions.hh"
#include <future>

//===================================================================
// Main program that constructs the grid and calls the driver routine
//===================================================================

int main(int argc, char** argv) {
	std::ofstream mylog;
	try {
		// Do some logging; it doesn't work yet.
		mylog.open("application.log");
		Dune::dinfo.attach(mylog);
		//Maybe initialize MPI
		Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
//		if (Dune::MPIHelper::isFake)
//			std::cout << "This is a sequential program." << std::endl;
//		else {
//			if (helper.rank() == 0)
//				std::cout << "parallel run on " << helper.size()
//						<< " process(es)" << std::endl;
//		}
        // Read the input file
		Params params;
		params.read_input(argc, argv);

		// sequential version
		if (helper.size() == 1) {
			const int dim = 2;

			// Construct 1D Bakhvalov grid which is able to resolve boundary layers
			// and use the Cartesian product of this 1D grid.
			std::vector<double> line;
			MGF<Params> mgf(params);  // mesh generating function
			mgf.double_side_interval(line);
			// Write down the grid points
			std::ofstream grid_pts("grid_pts.txt");
			for (auto & x : line)
				grid_pts << x << "\n";
			grid_pts.close();
			// make a tensor product grid
			Dune::array<std::vector<double>, dim> coords;
			for (auto& v : coords)
				v = line;
			typedef Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double, dim> > Grid;
			Grid grid(coords);

			grid.globalRefine(params.level);
			typedef Grid::LeafGridView GV;
			const auto& gv = grid.leafGridView();

			// Numerical models launch as parallel jobs.
			std::vector<std::future<int>> rets;  // execution times (in secs) for numerical models (parallel jobs)
			int anC = 0, anV=0, anN=0, an1=0; // execution times (in secs) for analytic model (serial)
			for(unsigned int i = 0; i < params.size; ++i){
				if(params.simulation[i] && i != Params::analytic_const
				                        && i != Params::analytic_new
				                        && i != Params::analytic_new1
						                && i != Params::analytic_var)
				{
				   params.model = static_cast<Params::Model>(i);
				   std::cout << "Model " << i << " launched.\n";
    			   rets.push_back(std::async(std::launch::async, driver<GV, Params>, gv, params));
				}
			}
			// Analytic jobs launch sequentially.
			for(unsigned int i = 0; i < params.size; ++i){
			    if(params.simulation[i] && (i == Params::analytic_const || i == Params::analytic_var
			    		                   || i == Params::analytic_new || i == Params::analytic_new1)){
				    params.model = static_cast<Params::Model>(i);
				    std::cout << "Solving " << params.simulation_names[i] << " sequentially.\n";
			    	int a = driver(gv, params);
			    	if(i == Params::analytic_const) anC = a;
			    	else if(i == Params::analytic_new) anN = a;
			    	else if(i == Params::analytic_new1) an1 = a;
			    	else anV = a;
			    }
			}
			// wait for all threads to complete before calling gnu_compare_c().
			std::cout << "Waiting...\n";
			for (unsigned int i = 0; i < rets.size(); ++i) {
				rets[i].wait();
			}
			// print execution times
			unsigned int ii = -1;
			for (unsigned int i = 0; i < params.size; ++i) {
				if (params.simulation[i] && i != Params::analytic_const
		                                 && i != Params::analytic_new
		                                 && i != Params::analytic_new1
		                                 && i != Params::analytic_var) {
					ii++;
					std::cout << " Time for " << params.simulation_names[i]
							<< " is = " << rets[ii].get() << " sec\n";
				}
			}
			if(params.simulation[Params::analytic_const])
				std::cout << " Time for " << params.simulation_names[Params::analytic_const]
											<< " is = " << anC << " sec\n";
			if(params.simulation[Params::analytic_var])
				std::cout << " Time for " << params.simulation_names[Params::analytic_var]
											<< " is = " << anV << " sec\n";
			if(params.simulation[Params::analytic_new])
				std::cout << " Time for " << params.simulation_names[Params::analytic_new]
											<< " is = " << anN << " sec\n";
			if(params.simulation[Params::analytic_new1])
				std::cout << " Time for " << params.simulation_names[Params::analytic_new1]
											<< " is = " << an1 << " sec\n";


			// Gnuplot control file for displaying the solution output is given only in 1D
			if(dim == 1) aux::gnu_output_solution(params);
			// Write gnuplot control file for displaying the fluxes.
			aux::gnu_compare_c(params);
		}

		Dune::dinfo.detach();
		mylog.close();
	} catch (Dune::Exception &e) {
		Dune::dinfo.detach();
		mylog.close();
		std::cerr << "Dune reported error: " << e << std::endl;
		return 1;
	} catch (std::exception &e) {
		Dune::dinfo.detach();
		mylog.close();
		std::cerr << "Std reported error: " << e.what() << std::endl;
		return 2;
	} catch (...) {
		Dune::dinfo.detach();
		mylog.close();
		std::cerr << "Unknown exception thrown!" << std::endl;
		return 3;
	}
}
