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
#include <vector>
#include <functional>

#include "imbibitionFunctions.hh"
#include "vangenuchten.hh"
//#include <dune/common/parametertree.hh>

#include <string>
#include <ctime>
#include <stdexcept>
#include <cstdlib>
#include <array>
#include <utility>
#include <string>

#include <boost/tokenizer.hpp>

namespace aux{

/// Date and time strig that is used as a name of the simulation folder.
std::string date_time(){
    std::time_t t = std::time(NULL);
    char mbstr[128];
    if (std::strftime(mbstr, sizeof(mbstr), "%F_%T", std::localtime(&t))) return mbstr;
    throw std::runtime_error("Date and time string error!");
}

/// Calculate data min and max for gnuplot command file
std::pair<double, double> min_max(std::string const & file_name, int colon){
   std::pair<double, double> tmp={1.0E100, -1.0E100}; // min, max

   std::ifstream in(file_name);
//   std::cout << file_name << std::endl;
   if(!in) throw std::runtime_error("Cannot open the file "+file_name);
   std::string line;
   while(std::getline(in, line)){
	   if(line[0] =='#') continue;
	   boost::char_separator<char> sep(" ");
	   boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
//	   auto size = tokens.end() - tokens.begin();
//	   if(colon > size) throw std::runtime_error("Too few collons in file "+file_name);
	   auto it = tokens.begin();
	   int i=0;
	   for(; it != tokens.end() && i < colon; ++i) ++it;
	   if(colon != i) throw std::runtime_error("Too few collons in file "+file_name);
//	   std::cout << *it << std::endl;
	   double value = std::stod(*it);
	   if(value < tmp.first) tmp.first = value;
	   if(value > tmp.second) tmp.second = value;
   }
   return tmp;
}

void create_dir(std::string const & dir_name) {
    std::string command = "mkdir -p "+dir_name;
    std::system(command.c_str());
}

template <typename Params>
void gnu_output_solution(Params const & params){
	if(!params.txtout) return;
	// count aktive simulations
	int total_cnt = 0;
	for(auto x : params.simulation) if(x) total_cnt++;

    const std::string & dir = params.date_and_time;
    std::string file = dir + "/solution.gnu";
    std::ofstream out(file);
    out << "set xrange [0:0.1]\n";
    out << "set yrange [0.:1.01]\n";
    out << "do for [ii=0:"<< params.Nout << "] {\n";
    int cnt = 0;
    for(unsigned int i = 0; i < params.size; ++i){
    	if(i == 0) out <<"   plot ";
    	if(params.simulation[i]){
    		cnt++;
    	    out << "  sprintf(\'"<< params.str_fname
			    << params.simulation_names[i] << "%d.txt\',ii) w l t \""
			    << params.simulation_names[i];
    		if(cnt != total_cnt)
    	       out << "\",\\\n";
    		else
    	       out << "\"\n";
    	}
    }
    out << "    pause 0.2\n";
    out << "}\n";
    out << "pause -1\n";
    out.close();
    std::cout << " To see the simulation results in gnuplot format run:\n  cd " << params.date_and_time
    		  << "; gnuplot solution.gnu\n";
}

template <typename Params>
void gnu_compare_c(Params const & params){
	int total_cnt = 0;
	for(auto x : params.simulation) if(x) total_cnt++;

	const std::string & dir = params.date_and_time;
	std::string file = dir + "/flux.gnu";
	std::ofstream out(file);
	out << "set xlabel \"time\"\n";
	out << "set ylabel \"Nonwetting source\"\n";
	out << "set key left center\n";
	out << "#set grid\n";
	out << "set xrange [0:1]\n";

	std::pair<double, double> tmp={1.0E100, -1.0E100}; // min, max
	for (unsigned int i = 0; i < params.size; ++i) {
		if (params.simulation[i]) {
	     	std::string name = params.str_sname  + params.simulation_names[i] + "-flux.txt";
	     	if(i == params.analytic_const){
	     	    auto tmp1 = min_max(name, 1);
	     	    if(tmp1.first < tmp.first) tmp.first = tmp1.first;
	     	    if(tmp1.second > tmp.second) tmp.second = tmp1.second;
	     	}
	     	else{
	     	    auto tmp1 = min_max(name, 1);
	     	    if(tmp1.first < tmp.first) tmp.first = tmp1.first;
	     	    if(tmp1.second > tmp.second) tmp.second = tmp1.second;
	     	    auto tmp2 = min_max(name, 2);
	     	    if(tmp2.first < tmp.first) tmp.first = tmp2.first;
	     	    if(tmp2.second > tmp.second) tmp.second = tmp2.second;

	     	}
		}
	}
	double dy = tmp.second - tmp.first;
	assert(dy >= 0.0);
	tmp.first -= 0.05*dy;
	tmp.second += 0.05*dy;
	out << "set yrange [" << tmp.first <<":" << tmp.second << "]\n";
	out << "#set terminal postscript eps color solid lw 3\n";
	out << "#set output \"clin_Q.eps\"\n";
	int cnt = 0;
	for (unsigned int i = 0; i < params.size; ++i) {
     	std::string name = params.simulation_names[i];
		if (i == 0)
			out << "   plot ";
		if (params.simulation[i]) {
			cnt++;
			if(i == params.analytic_const)
     			out << "\"" << name << "-flux.txt\" u 1:2 w l t \"" << name;
			else if(i == params.analytic_var)
     			out << "\"" << name << "-flux.txt\" u 1:2 w l t \"" << name;
			else{
			   out << "\"" << name << "-flux.txt\" u 1:2 w l t \"" << name + " vol int\",\\\n";
			   out << "\"" <<  name << "-flux.txt\" u 1:3 w l t \"" << name + " bdr int";
			}
			if (cnt != total_cnt)
				out << "\",\\\n";
			else
				out << "\"\n";
		}
	}
	out << "pause -1\n";
	out.close();
	std::cout
			<< " To see the flux comparison in constant linear case in gnuplot format run:\n  cd "
			<< params.date_and_time << "; gnuplot flux.gnu\n";
}




} // end of namespace aux




/** \brief Structure that gives all data necessary for simulation.
 *
 * This class reads the parameters from the input file. All data necessary for
 * a simulation are present in this class.
 *
 * *Data*
 *
 * -   type of model
 * -   functions \f$\alpha(S)\f$ and \f$\beta(S)\f$ and boundary condition function
 * -   permeability, porosity, \f$\delta\f$, domain length
 * -   grid generation parameters
 * -   time stepping parameters
 * -   output file parameters
 *
 * *Flux functions*
 * - FluxFunction = 0 selects \f$ \alpha(S) = aS(1-S)\f$. Parameter \f$a\f$ must be given.
 * - FluxFunction = 1 selects \f$ \alpha(S) = -\lambda_w(S) \lambda_n(S) p_c'(S)/\lambda(S) \f$
 *    for van Genuchten functions with given parameters.
 *
 * *Boundary function*
 * - The Dirichlet boundary condition can be selected by *BoundaryFunction* index. The functions
 *  are hard coded as lambdas.
 *
 */
struct Params{
	/// Enum constants describing different imbibition models.
	enum Model{
	  new_nonlinear=0, nonlinear, constant_linear, variable_linear, analytic_const, analytic_var, size
	};
 /**
  * Read all parameters from an input file. It must be called explicitly
  * since it is not called in the constructor.
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
    str_fname =""; //  input_data.get<std::string> ("Output.SimulationBaseName");  -- not needed
    Nout = tend/dtout;
    std::string str_model =  input_data.get<std::string> ("Model");
    set_simulation(str_model);
    function_index  = input_data.get<int>("BoundaryFunction");
    flux_funct_index = input_data.get<int>("FluxFunction");

    if(function_index < 0 or function_index >= ptfun.size())
       throw std::runtime_error("Wrong boundary function index! index = " + std::to_string(function_index));
    if(flux_funct_index < 0 or flux_funct_index > 1)
       throw std::runtime_error("Wrong flux function index! index = " + std::to_string(flux_funct_index));

    // read parameters for selected function
    if(flux_funct_index == 0){
    	// Need only parameter a.
        a  =  input_data.get<double>("AlphaFunction.Amplitude");
        aImbFun.init(a);  // this object will be used
        mean_alpha = aImbFun.beta(1.0);
    }
    else
    {
       double vgAlpha  =  input_data.get<double>      ("VanGenuchten.Alpha");
       double vgN      =  input_data.get<double>      ("VanGenuchten.N");
       Dumux::VanGenuchtenParams vgParams(vgAlpha, vgN);

       double muw =  input_data.get<double>("Fluids.WettingViscosity");
       double mun =  input_data.get<double>("Fluids.NonWettingViscosity");
       // this object will be used
       vgImbFun.init(vgParams, muw, mun);
       mean_alpha = vgImbFun.beta(1.0);
    }
    amin = mean_alpha * acom;
    // All simulation output goes to the folder named after current date and time.
    date_and_time = aux::date_time();
    // create folder if it does not exist.
    aux::create_dir(date_and_time);
    str_sname = "./" + date_and_time + "/" + str_fname;
    // copy the input file into the simulation folder for documentation purposes
    std::string command = "cp " + filename + " " + date_and_time +"/";
    std::system(command.c_str());
    return;
  }

   /// Constructor.
   Params(std::string const & file_name = "imbibition.input") : default_file_name(file_name)
   {
         // mean_alpha = integrate_alpha();
         // std::cout << "mean alpha  " << mean_alpha << std::endl;
          ptfun.push_back([](double t) { return 0.35 + 0.2* std::sin(2*M_PI*t); });
          ptfun.push_back([](double t) { return ((0.9-0.1)/M_PI)*std::atan((2*t-1)/2.5)+(0.9+0.1)/2; });
          ptfun.push_back([](double t) { return -((0.9-0.1)/M_PI)*std::atan((2*t-1)/0.1)+(0.9+0.1)/2; });
          ptfun.push_back([](double t) { return 0.35 - 0.2* std::sin(2*M_PI*t); });
          ptfun.push_back([](double t) { return std::max(std::min(0.5 + 0.51* std::sin(2*M_PI*t), 1.0), 0.0); });
          ptfun.push_back([](double t) { return 1.0; });

          for(unsigned int i=0; i < size; ++i) simulation[i] = false;
          simulation_names[new_nonlinear] = "n_nlin";
          simulation_names[nonlinear] = "nlin";
          simulation_names[constant_linear] = "clin";
          simulation_names[variable_linear] = "vlin";
          simulation_names[analytic_const] = "anac";
          simulation_names[analytic_var] = "anav";
   }

   Dune::ParameterTree input_data;
   // Type of model to solve. It mus the set before calling the driver.
   Model model = Model::size;
   ///  \f$\alpha(S)\f$  nonlinear diffusivity coefficient.
   double alpha(double u) const
   {
     if(flux_funct_index == 0) return aImbFun.alpha(u);
     return vgImbFun.alpha(u);
   }

   ///  \f$\beta(S) = \int_0^S \alpha(u) du\f$ .
   double beta(double u) const
   {
     if(flux_funct_index == 0) return aImbFun.beta(u);
     return vgImbFun.beta(u);
   }

   /// \f$\alpha(S)\f$  nonlinear diffusivity coefficient that is cut-off
   /// in order to remain strictly positive. Probably not needed.
   double alpha_reg(double u) const {
     double a = alpha(u);
     if(a < amin) a = amin;
     return a;
   }
   /// Saturation boundary condition on the matrix block boundary
   double bdry(double t) const { return ptfun[function_index](t);  }
   /// return boundary function (needed for analytic solution)
   std::function<double(double)> bdry_fun() const {return ptfun[function_index]; }

   // Constants
   double a = 0.0;       ///< amplitude of the artificial alpha function
   double acom = 0.0;    ///< alpha cut off multiplier
   double mean_alpha = 0.0;  ///<   \f$\int_0^1 \alpha(s) ds\f$
   // Porous media
   double k = 0.0;       ///< permeability
   double poro = 0.0;    ///< porosity
   // Grid generation parameters
   double delta = 0.0;   ///< \f$\delta\f$ parameter
   double q = 0.0;       ///< Bakhvalov grid generation parameter.
   double sigma = 0.0;   ///< Bakhvalov grid generation parameter.
   // grid
   double L = 0.0; ///<  domain side length
   int    N = 0; ///< number of points per side
   int    level = 0; ///< grid refinement level
   // time stepping
   double dt = 0.0;    ///< \f$\Delta t\f$
   double dtmax = 0.0; ///< maximum \f$\Delta t\f$
   double dtout = 0.0; ///< time step for output operation
   double tend = 0.0;  ///< final time of the simulation
   // output
   bool vtkout = false; ///< Do output VTK files
   bool txtout = false; ///< Do output GNUPLOT TXT files
   std::string str_fname; ///< simulation base name without the directory part
   std::string str_sname; ///< simulation base name with the directory part
   std::string date_and_time; ///< folder name for the output
   std::array<bool,size> simulation;  ///< simulations to make (true/false flags)
   std::array<std::string,size> simulation_names; ///< simulation names corresponding to simulations
   int Nout = 0; ///< Output files will be numbered from 0 to Nout (inclusive)
private:
   std::string default_file_name;
   std::vector<std::function<double(double)>> ptfun;
   unsigned int function_index = 0;
   unsigned int flux_funct_index = 0;
   double amin = 0.0; //alpha(0.5)/20;
   // implementations of alpha-functions
   ArtifImbibitionFunctions aImbFun;
   RealImbibitionFunctions<Dumux::VanGenuchten> vgImbFun;

   void set_simulation(std::string const & sim){
	   for(auto x : sim){
		   switch(x){
		   case 'c':
		   case 'C': simulation[constant_linear] = true;
		             break;
		   case 'v':
		   case 'V': simulation[variable_linear] = true;
		             break;
		   case 'n':
		   case 'N': simulation[nonlinear] = true;
		             break;
		   case 'm':
		   case 'M': simulation[new_nonlinear] = true;
		             break;
		   case 'a':
		   case 'A': simulation[analytic_const] = true;
		             break;
		   case 'b':
		   case 'B': simulation[analytic_var] = true;
		             break;

		   }
	   }
   }
};



#endif /* SRC_PARAMETERS_HH_ */
