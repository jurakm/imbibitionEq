// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

 \mainpage Solution of nonlinear imbibition equation and two of its linearizations.

 \section Problem

   We solve the imbibition equation and its two possible linearizations. 
 * We calculate  the solution and the matrix-fracture source term.
   Imbibition equation is a nonlinear parabolic equation:
   \f[
       \Phi_m \frac{\partial S}{\partial t} - \delta^2 k_m \Delta \beta_m(S) = 0 \quad {\rm in }\;  Y,\; t>0
   \f]
   \f[  S = {\cal P}(S^f)\quad{\rm on } \;\partial Y,\; t > 0 \f]
   \f[  S = S_0 \quad {\rm on}\; Y,\; t = 0 \f]
 * Description:
 * - \f$\Phi_m\f$ is the matrix block porosity. 
 * - \f$\delta\f$ is relative fracture thickness. That is, if the matrix block size is \f$l\f$ (not present here)
 *  then the fracture thickness is \f$l\delta\f$. 
 * - \f$S\f$ is the wetting phase saturation.
 * - \f$Y = (0,1)^d\f$. Actually, it should be  \f$(\delta/2,1-\delta/2)^d\f$ but we have neglected this for now.
 * 
    We also solve two linearized versions of the imbibition equation and
 * associated matrix-fracture flux (which is our main interest). There are several ways to
 * linearize the nonlinear imbibition problem:
 *  \f[
   \Phi_m \frac{\partial S}{\partial t} - \delta^2 k_m \Delta \beta_m(S) = 0 \quad {\rm in }\;  Y,\; t>0
   \f]
   \f[  S = g(t)\quad{\rm on } \;\partial Y,\; t > 0 \f]
   \f[  S = g(0)\quad {\rm on}\; Y,\; t = 0, \f]
 * where \f$Y = (0,1)^d\f$.
 *
 * Note. We assume fixed domain \f$Y\f$ (independent of \f$\delta\f$)
 * since \f$\delta\f$ is small and we can replace  \f$Y^\delta = (\delta/2,1-\delta/2)^d\f$
 * by \f$(0,1)^d\f$.
 *
 * All approximations are of the form
 * \f[
       \Phi_m \frac{\partial S}{\partial t} - \delta^2 k_m a(t)\Delta S = 0 \quad {\rm in }\;  Y,\; t>0
   \f]
   \f[  S = g(t)\quad{\rm on } \;\partial Y,\; t > 0 \f]
   \f[  S = g(0)\quad {\rm on}\; Y,\; t = 0, \f]
 * where for \f$a(t)\f$ we have following possibilities:
 *       - Constant approximation [anac]: \f$ a(t) = \bar\alpha_m,\; \tau(t) = t. \f$
 *       - Variable approximation [anav]: \f$ a(t) = \alpha_m(g(t)),\; \tau(t) = \int_0^t a(u)du. \f$
 *       - Variable approximation [ana_n]:
 *   \f$ a(t) = (\beta_m(g(t))-\beta_m(g(t) -\theta (g(t) -g(0)))/[\theta (g(t)-g(0))],\;
 *        \tau(t) = \int_0^t a(u)du. \f$
 *
 * where \f$\alpha_m(S) = \beta_m'(S)\f$ and  \f$\bar\alpha_m = \int_0^1 \alpha_m(s)\, ds\f$.
 * The function \f$a(t) \f$ is eliminated by passing to a new time variable
 * \f$\tau(t) = \int_0^t a(u)du\f$ (except in the first case where \f$a(t) =1\f$).

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
#include <future>

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
#include <dune/pdelab/backend/istlvectorbackend.hh>
//#include <dune/pdelab/backend/istl.hh> -- nova verzija pdelaba 
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include "parameters.hh"
//#include "read_input.hh"
#include "driver.hh"
#include "mgf.hh"
#include "linear_analytic.hh"
#include "imbibitionFunctions.hh"


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
		if (helper.size() != 1)
            throw std::runtime_error("Please launch program in a sequential mode. Aborting.");
            
        const int dim = 1;

        // Construct 1D Bakhvalov grid which is able to resolve boundary layers
        // and use the Cartesian product of this 1D grid.
        std::vector<double> line;
        MGF<Params> mgf(params);  // mesh generating function
        mgf.double_side_interval(line);
        // Write down the grid points
        std::ofstream grid_pts(params.date_and_time+"/grid_pts.txt");
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

        // Set O2SCL error handler 
        auto tmp = o2scl::err_hnd; 
        auto new_err_hnd = new o2scl::err_hnd_cpp();
            
        // Numerical models launch as parallel jobs.
        std::vector<std::future<int>> rets;  // execution times (in secs) for numerical models (parallel jobs)
//		int anC = 0, anV=0, anN=0, an1=0; // execution times (in secs) for analytic model (serial)
        for(unsigned int i = 0; i < params.size; ++i){
            if(params.simulation[i])
            {
                params.model = static_cast<Params::Model>(i);
                std::cout << "Model " << params.simulation_names[i] << " launched.\n";
                rets.push_back(std::async(std::launch::async, driver<GV, Params>, gv, params));
                
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
            if (params.simulation[i])
            {
                ii++;
                std::cout << " Time for " << params.simulation_names[i]
                        << " is = " << rets[ii].get() << " sec\n";
            }
        }

        // Gnuplot control file for displaying the solution output is given only in 1D
        if(dim == 1) aux::gnu_output_solution(params);
        // Write gnuplot control file for displaying the fluxes.
        aux::gnu_compare_c(params);
            
        // O2SCL. Restore old error handler.
        o2scl::err_hnd = tmp;  
        delete new_err_hnd;
            
		
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
