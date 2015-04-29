/*
 * parameters.hh
 *
 *  Created on: 22. sij 2015.
 *      Author: jurak
 */

#ifndef SRC_PARAMETERS_HH_
#define SRC_PARAMETERS_HH_

#include <cmath>
#include <cstdlib>

#include "imbibitionFunctions.hh"
#include "vangenuchten.hh"
//#include <dune/common/parametertree.hh>

enum class Model{
  nonlinear, constant_linear, variable_linear, all, size
};

// Boundary value functions
double bdry_0(double t) { return 0.35 + 0.2* std::sin(2*M_PI*t); }
double bdry_1(double t) { return ((0.9-0.1)/M_PI)*std::atan((2*t-1)/2.5)+(0.9+0.1)/2; }
double bdry_2(double t) { return -((0.9-0.1)/M_PI)*std::atan((2*t-1)/0.1)+(0.9+0.1)/2; }
double bdry_3(double t) { return 0.35 - 0.2* std::sin(2*M_PI*t); }
double bdry_4(double t) { return 0.5 + 0.49* std::sin(2*M_PI*t); }
double bdry_5(double t) { return 1.0; }

/** Structure that gives all data necessary for simulation:
 */
template <typename ParameterTree>
struct Params{

 /**
  * Read all parameters from an input file.
  */
  void read_input (int argc, char** argv)
  {
    std::string filename = default_file_name;

    if (argc > 1) filename = argv[1];

    try {
        Dune::ParameterTreeParser::readINITree (filename, input_data);
      }
    catch (...) {
        std::cerr << "The configuration file \"" << filename << "\" "
  	  "could not be read. Exiting..." << std::endl;
        std::exit(1);
      }

    k         =  input_data.get<double>      ("Permeability");
    poro      =  input_data.get<double>      ("Porosity");
    delta     =  input_data.get<double>      ("Delta");
    acom      =  input_data.get<double>      ("AlphaCutOffMultiplier");
    level     =  input_data.get<int>         ("Refinement.Level");
    N         =  input_data.get<int>         ("Grid.NPoints");
    q         =  input_data.get<double>      ("Grid.Q");
    L         =  input_data.get<double>      ("Grid.Length");
    sigma     =  input_data.get<double>      ("Grid.Sigma");
    dt        =  input_data.get<double>      ("Time.Dt");
    dtmax     =  input_data.get<double>      ("Time.DtMax");
    dtout     =  input_data.get<double>      ("Time.DtOut");
    tend      =  input_data.get<double>      ("Time.Final");
    vtkout    =  input_data.get<int>         ("Output.VTK");
    txtout    =  input_data.get<int>         ("Output.TXT");
    str_sname =  input_data.get<std::string> ("Output.SimulationBaseName");

    std::string str_model =  input_data.get<std::string> ("Model");
    if(str_model == std::string("constant_linear"))
      model = Model::constant_linear;
    else if(str_model == std::string("variable_linear"))
      model = Model::variable_linear;
    else if(str_model == std::string("nonlinear"))
      model = Model::nonlinear;
    else if(str_model == std::string("all"))
      model = Model::all;
    else throw std::runtime_error("Model from input file is unknown");

    a  =  input_data.get<double>("AlphaFunction.Amplitude");
    function_index  = input_data.get<int>("BoundaryFunction");
    flux_funct_index = input_data.get<int>("FluxFunction");

    if(function_index < 0 or function_index > 5)
       throw std::runtime_error("Wrong boundary function index! index = " + std::to_string(function_index));
    if(flux_funct_index < 0 or flux_funct_index > 1)
       throw std::runtime_error("Wrong flux function index! index = " + std::to_string(flux_funct_index));


    double vgAlpha  =  input_data.get<double>      ("VanGenuchten.Alpha");
    double vgN      =  input_data.get<double>      ("VanGenuchten.N");
    Dumux::VanGenuchtenParams vgParams(vgAlpha, vgN);

    double muw =  input_data.get<double>("Fluids.WettingViscosity");
    double mun =  input_data.get<double>("Fluids.NonWettingViscosity");
    vgImbFun.init(vgParams, muw, mun);
    aImbFun.init(a);

    if(flux_funct_index == 0)
      mean_alpha = aImbFun.beta(1.0);
    else
      mean_alpha = vgImbFun.beta(1.0);

    amin = mean_alpha * acom;
  }

   Params(std::string const & file_name = "imbibition.input") : default_file_name(file_name)
   {
         // mean_alpha = integrate_alpha();
         // std::cout << "mean alpha  " << mean_alpha << std::endl;
          ptfun[0] = bdry_0;
          ptfun[1] = bdry_1;
          ptfun[2] = bdry_2;
          ptfun[3] = bdry_3;
          ptfun[4] = bdry_4;
          ptfun[5] = bdry_5;
   }

   Dune::ParameterTree input_data;
   // Type of model to solve
   Model model = Model::size;
   // Nonlinearity
   double alpha(double u) const
   {
     if(flux_funct_index == 0) return aImbFun.alpha(u);
     return vgImbFun.alpha(u);
   }

   double beta(double u) const
   {
     if(flux_funct_index == 0) return aImbFun.beta(u);
     return vgImbFun.beta(u);
   }

   double alpha_reg(double u) const {
     double a = alpha(u);
     if(a < amin) a = amin;
     return a;
   }
   // Boundary condition
   double bdry(double t) const { return ptfun[function_index](t);  }


   // Constants
   double a = 0.0;       // amplitude of the artificial alpha function
   double acom = 0.0;    // alpha cut off multiplier
   double mean_alpha = 0.0;
   // Porous media
   double k = 0.0;       // permeability
   double poro = 0.0;    // porosity
   // Grid generation parameters
   double delta = 0.0;   // small parameter
   double q = 0.0;
   double sigma = 0.0;
   // grid
   double L = 0.0; // side size
   int    N = 0; // number of points per size
   int    level = 0; // refinement level
   // time stepping
   double dt = 0.0;
   double dtmax = 0.0;
   double dtout = 0.0;
   double tend = 0.0;
   // output
   bool vtkout = false;
   bool txtout = false;
   std::string str_sname; // simulation name
private:
   std::string default_file_name;
   std::array<double(*)(double),6> ptfun;
   unsigned int function_index = 0;
   unsigned int flux_funct_index = 0;
   double amin = 0.0; //alpha(0.5)/20;

   ArtifImbibitionFunctions aImbFun;
   RealImbibitionFunctions<Dumux::VanGenuchten> vgImbFun;
};



#endif /* SRC_PARAMETERS_HH_ */
