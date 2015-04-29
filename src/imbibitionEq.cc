// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    
    \brief Parabolička jednadžba diskretizirana konformnom metodom KE u prostoru 
          i implicitnom Eulerovom metodom u vremenu. 
    */
#ifdef HAVE_CONFIG_H
#include "config.h"     
#endif

#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/timer.hh>
// Za čitanje .input datoteke
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>


#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>
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
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/instationary/onestep.hh>

#include "parameters.hh"
//#include "read_input.hh"
#include "driver.hh"
#include "mgf.hh"

#include <future>

//===================================================================
// Main program that constructs the grid and calls the driver routine
//===================================================================

/**
 * Solution of nonlinear imbibition equation and two of its linearizations.
 */
int
main (int argc, char** argv)
{
  std::ofstream mylog;
  try
    {
      mylog.open("application.log");
      Dune::dinfo.attach(mylog);
      //Maybe initialize MPI
      Dune::MPIHelper& helper = Dune::MPIHelper::instance (argc, argv);
      if (Dune::MPIHelper::isFake)
	std::cout << "This is a sequential program." << std::endl;
      else
	{
	  if (helper.rank () == 0)
	    std::cout << "parallel run on " << helper.size () << " process(es)"
		<< std::endl;
	}

      using Parameters = Params<Dune::ParameterTree>;
      Parameters params;
      params.read_input(argc, argv);


      bool test = false;                 // test integration or not
      if (argc > 4)
	test = true;

      // sequential version
      if (helper.size () == 1)
	{
	  const int dim = 1;

	  // We construct 1D Backhvalov grid which is able to resolve boundary layers
	  // and then we use the Cartesian product of this 1D grid.
	  std::vector<double> line;
	  MGF<Parameters> mgf(params);  // mesh generating function
	  mgf.double_side_interval(line);
	  // Write down the grid points
	  std::ofstream grid_pts("grid_pts.txt");
	  for(auto & x : line) grid_pts << x << "\n";
	  grid_pts.close();
          // make a tensor product grid
	  Dune::array<std::vector<double>, dim> coords;
	  for (auto& v : coords)
	    v = line;
	  typedef Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double, dim> > Grid;
	  Grid grid (coords);

	  grid.globalRefine (params.level);
          typedef Grid::LeafGridView GV;
	  const auto& gv = grid.leafGridView ();

	  if (params.model == Model::all)
	    {
	      params.model = Model::constant_linear;
	      std::async (std::launch::async, driver<GV, Parameters>, gv, params);
	      params.model = Model::variable_linear;
	      std::async (std::launch::async, driver<GV, Parameters>, gv, params);
	      params.model = Model::nonlinear;
	      std::async (std::launch::async, driver<GV, Parameters>, gv, params);
	    }
	  else
	    driver (gv, params);
	  // Test the routine for volume integration
	}

     Dune::dinfo.detach();
     mylog.close();
    }
  catch (Dune::Exception &e)
    {
      Dune::dinfo.detach();
      mylog.close();
      std::cerr << "Dune reported error: " << e << std::endl;
      return 1;
    }
  catch (std::exception &e)
     {
     Dune::dinfo.detach();
     mylog.close();
       std::cerr << "Std reported error: " << e.what() << std::endl;
       return 2;
     }
  catch (...)
    {
      Dune::dinfo.detach();
      mylog.close();
      std::cerr << "Unknown exception thrown!" << std::endl;
      return 3;
    }
}
