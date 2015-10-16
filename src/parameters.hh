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

#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

namespace aux{

std::string sec_to_string(int sec){
	std::string res;
	if(sec < 60) res = std::to_string(sec);
	else if(sec < 60*60){
		int s = sec % 60;
		int m = sec/60;
		res = std::to_string(m) + " min "+ std::to_string(s) + " sec";
	}
	else{
		int s = sec % 60;
		int m = sec/60;
		m = m % 60;
		int h = m/60;
		res = std::to_string(h) + " hours " + std::to_string(m) + " min "+ std::to_string(s) + " sec";
	}
	return res;
}

/// Date and time string that is used as a name of the simulation folder.
std::string date_time(){
    std::time_t t = std::time(NULL);
    char mbstr[128];
    if (std::strftime(mbstr, sizeof(mbstr), "%b-%d.%H.%M.%S", std::localtime(&t))) return mbstr;
    throw std::runtime_error("Date and time string error!");
}

/** Calculate data minimum and maximum of a given column in a file.
 *  The file consists of the columns of data and possible comment lines
 *  that start with the sign '#' in the first column. Comment lines are ignored.
 *  @param file_name = file name
 *  @param colon = colon index, starting from zero.
 */
std::pair<double, double> min_max(std::string const & file_name, int colon)
		{
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
//	   if(colon > size) throw std::runtime_error("Too few colons in file "+file_name);
	   auto it = tokens.begin();
	   int i=0;
	   for(; it != tokens.end() && i < colon; ++i) ++it;
	   if(colon != i) throw std::runtime_error("Too few colons in file "+file_name);
//	   std::cout << *it << std::endl;
	   double value = std::stod(*it);
	   if(value < tmp.first) tmp.first = value;
	   if(value > tmp.second) tmp.second = value;
   }
   return tmp;
}

/** Create directory (linux specific). */
void create_dir(std::string const & dir_name) {
    std::string command = "mkdir -p "+dir_name;
    std::system(command.c_str());
}

/// Write gnuplot command file to see an animation of the solution(s).
template <typename Params>
void gnu_output_solution(Params const & params){
	if(!params.txtout) return;
	// count active simulations
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


/// Write gnuplot command file to see the fluxes. Plot certain combination
/// of fluxes.
// Params::new_nonlinear,   Params::nonlinear,    Params::constant_linear,
// Params::variable_linear, Params::variable_new, Params::analytic_const,
// Params::analytic_var,    Params::analytic_new, Params::analytic_new1
template <typename Params>
void gnu_compare_c(Params const & params){
	params.gnu_compare_c();
	params.gnu_compare_c({Params::new_nonlinear, Params::nonlinear}, "-nlin");
	params.gnu_compare_c({Params::nonlinear,
							  Params::constant_linear,
							  Params::variable_linear,
							  Params::variable_new}, "-cmp");
	params.gnu_compare_c({Params::constant_linear,Params::analytic_const}, "-const");
	params.gnu_compare_c({Params::analytic_var, Params::analytic_new,
		Params::variable_linear, Params::variable_new, Params::analytic_new1}, "-var");

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
	/// Constants describing different imbibition models.
	enum Model{
	  new_nonlinear=0, nonlinear, constant_linear, variable_linear, variable_new,
	  analytic_const, analytic_var, analytic_new, analytic_new1, size
	};
 /**
  * Read all parameters from an input file. It must be called explicitly
  * since it is not called in the constructor.
  */
  void read_input (int argc, char** argv);


  /// Constructor. Defines boundary functions (evolution of fracture saturations)
  /// and simulation names.
   explicit	Params(std::string const & file_name = "imbibition.input");

   Dune::ParameterTree input_data;
   /// Type of model to solve. It must the set before calling the driver.
   Model model = Model::size;
   ///  \f$\alpha_m(S)\f$  nonlinear diffusivity coefficient in the matrix.
   double alpha(double u) const
   {
     if(flux_funct_index == 0) return aImbFun.alpha(u);
     return vgImbFunMatrix.alpha(u);
   }

   ///  \f$\beta_m(S) = \int_0^S \alpha_m(u) du\f$ .
   double beta(double u) const
   {
     if(flux_funct_index == 0) return aImbFun.beta(u);
     return vgImbFunMatrix.beta(u);
   }

   /// Matrix capillary pressure function.
   double pc_matrix(double sw) const {
	  if(flux_funct_index == 0) return sw;
	  return vgImbFunMatrix.pc( sw );
   }

   /// Fracture capillary pressure function.
   double pc_fracture (double sw) const {
   	  if(flux_funct_index == 0) return sw;
   	  return vgImbFunFracture.pc( sw );
   }

   /// Boundary transfer function transforming fracture saturation to matrix saturation.
   double bdry_transfer(double s) const{
	   if(flux_funct_index == 0) return s;
	   return vgImbFunMatrix.sw( vgImbFunFracture.pc( s ) );
   }

	/// Saturation boundary condition on the matrix block boundary.
	double bdry(double t) const {
		if (flux_funct_index == 0)
			return ptfun[function_index](t);
		return bdry_transfer(ptfun[function_index](t));
	}

	/** Linearized diffusivity coefficient. For different models we have
	 *  different functions.
	 *  */
	double a_g(double t) const {
		static const double TOL = 1.0e-5;
		double val = 1.0; // good value for Params::new_nonlinear -- important.
		if (model == analytic_const or model == constant_linear)
			val = mean_alpha;
		else if (model == analytic_var or model == variable_linear)
			val = alpha(bdry(t));
		else if (model == analytic_new or model == variable_new) {
			// we have the best results with theta = 0.75.
			const double Yt = bdry(t);
			const double dS = Yt - bdry(0.0);
			if (std::abs(dS * theta) > TOL) {
				val = (beta(Yt) - beta(Yt - dS * theta)) / (dS * theta);
			} else
				val = alpha(Yt);
		} else if (model == Params::analytic_new1) {
			const double Yt = bdry(t);
			const double Yt0 = (t > dt_bdry) ? bdry(t - dt_bdry) : bdry(0.0);
			const double dS = Yt - Yt0;
			if (std::abs(dS) > TOL) {
				val = (beta(Yt) - beta(Yt0)) / dS;
			} else
				val = alpha(Yt);
		}
		return val;
	}
//	const ImbibitionFunctions * const imbib_fun() const {
//		if (flux_funct_index == 0)
//			return &aImbFun;
//		return &vgImbFunMatrix;
//	}
	// Constants
	double a = 0.0;       ///< amplitude of the artificial alpha function
   double mean_alpha = 0.0;  ///<   \f$\int_0^1 \alpha(s) ds\f$
   // Porous media
   double k = 0.0;       ///< permeability
   double poro = 0.0;    ///< porosity
   double theta = 1.0;   ///< a factor in calculating the mean value of diffusion coefficient
   double dt_bdry = 0.0; ///< a dt in calculating the mean value of diffusion coefficient
   // Grid generation parameters
   double delta = 0.0;   ///< \f$\delta\f$ parameter
   double scaled_delta = 0.0; ///< \f$\delta\sqrt{k\alpha_m/\Phi}\f$
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


   /**
      * Write gnuplot command file to see the fluxes.
      * @param show_sim = set of fluxes to show is they were simulated. By default show
      *                   all simulated fluxes.
      * @param add_to_name = string to add to a base file name to distinguish different
      *                      file names. Default (for all simulations) is empty string.
      */
     void gnu_compare_c(std::set<int> const & show_sim =
          {new_nonlinear, nonlinear, constant_linear, variable_linear, variable_new,
         		 analytic_const, analytic_var, analytic_new, analytic_new1},
  			 std::string const & add_to_name ="") const;

private:
   std::string default_file_name;
   std::vector<std::function<double(double)>> ptfun;
   unsigned int function_index = 0;
   unsigned int flux_funct_index = 0;
   // implementations of alpha-functions
   ArtifImbibitionFunctions aImbFun;
   RealImbibitionFunctions<Dumux::VanGenuchten> vgImbFunMatrix;
   RealImbibitionFunctions<Dumux::VanGenuchten> vgImbFunFracture;

   void set_simulation(std::string const & sim){
	   for(auto x : sim){
		   switch(x){
		   case 'c':
		   case 'C': simulation[constant_linear] = true;
		             break;
		   case 'v':
		   case 'V': simulation[variable_linear] = true;
		             break;
		   case 'z':
		   case 'Z': simulation[variable_new] = true;
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
		   case 'e':
		   case 'E': simulation[analytic_new] = true;
		             break;
		   case 'f':
		   case 'F': simulation[analytic_new1] = true;
		             break;
		   }
	   }
   }
   /** Plot all given functions of the model.
    * @param n = number of points in tables. If not given a default is used.
    * */
   void plot_functions(unsigned int n = 0) const;

};

Params::Params(std::string const & file_name) :	default_file_name(file_name) {
		ptfun.push_back([](double t) {return 0.35 + 0.2* std::sin(2*M_PI*t);});
		ptfun.push_back(
				[](double t) {return ((0.9-0.1)/M_PI)*std::atan((2*t-1)/2.5)+(0.9+0.1)/2;});
		ptfun.push_back(
				[](double t) {return -((0.9-0.1)/M_PI)*std::atan((2*t-1)/0.1)+(0.9+0.1)/2;});
		ptfun.push_back([](double t) {return 0.35 - 0.2* std::sin(2*M_PI*t);});
		ptfun.push_back(
				[](double t) {return std::max(std::min(0.5 + 0.51* std::sin(2*M_PI*t), 1.0), 0.0);});
		ptfun.push_back([](double t) {return 0.05 + std::min(t,0.9);});
		ptfun.push_back([](double t) {return 0.05 + std::min(t/10.0, 0.9);});

		for (unsigned int i = 0; i < size; ++i)
			simulation[i] = false;
		simulation_names[new_nonlinear] = "n_nlin";
		simulation_names[nonlinear] = "nlin";
		simulation_names[constant_linear] = "clin";
		simulation_names[variable_linear] = "vlin";
		simulation_names[variable_new] = "vlin_n";
		simulation_names[analytic_const] = "anac";
		simulation_names[analytic_var] = "anav";
		simulation_names[analytic_new] = "ana_n";
		simulation_names[analytic_new1] = "ana_1";
	}


/**
 * Read all parameters from an input file. It must be called explicitly
 * since it is not called in the constructor.
 */
 void Params::read_input (int argc, char** argv)
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
   theta     =  input_data.get<double>      ("Theta");
   dt_bdry   =  input_data.get<double>      ("DtBdry");
//    acom      =  input_data.get<double>      ("AlphaCutOffMultiplier");
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
      double vgAlphaMat  =  input_data.get<double>      ("Matrix-VanGenuchten.Alpha");
      double vgNMat   =  input_data.get<double>      ("Matrix-VanGenuchten.N");
      Dumux::VanGenuchtenParams vgParamsMatrix(vgAlphaMat, vgNMat);
      double vgAlphaFr  =  input_data.get<double>      ("Fracture-VanGenuchten.Alpha");
      double vgNFr      =  input_data.get<double>      ("Fracture-VanGenuchten.N");
      Dumux::VanGenuchtenParams vgParamsFracture(vgAlphaFr, vgNFr);

      double muw =  input_data.get<double>("Fluids.WettingViscosity");
      double mun =  input_data.get<double>("Fluids.NonWettingViscosity");
      // this object will be used
      vgImbFunMatrix.init(vgParamsMatrix, muw, mun);
      mean_alpha = vgImbFunMatrix.beta(1.0);
      vgImbFunFracture.init(vgParamsFracture, muw, mun);
   }
//    amin = mean_alpha * acom;
   // All simulation output goes to the folder named after current date and time.
   date_and_time = aux::date_time();
   // create folder if it does not exist.
   aux::create_dir(date_and_time);
   str_sname = "./" + date_and_time + "/" + str_fname;
   // copy the input file into the simulation folder for documentation purposes
   std::string command = "cp " + filename + " " + date_and_time +"/";
   std::system(command.c_str());

   // plot alpha and beta functions for gnuplot
	plot_functions();
	scaled_delta = delta * std::sqrt(k*mean_alpha/poro);
	std::cout << "Scaled delta = " << scaled_delta << std::endl;
   return;
 }

void Params::plot_functions(unsigned int n) const
{
	const ImbibitionFunctions * pfun = nullptr;
	if (flux_funct_index == 0) pfun = &aImbFun;
	else                       pfun = &vgImbFunMatrix;

    if(n == 0) n = pfun->NofPts();
    // functions are defined on [0,1]
    double hh = 1.0/n;

    const std::string & dir = date_and_time;
    std::string file_name = "functions.txt";
    std::string full_file_name = dir+"/functions.txt";

    std::ofstream file (full_file_name);
    file << "#    S_w        alpha(S_w)        beta(S_w)     bdry_trans(S_w)\n";
    for (unsigned int i = 0; i <= n; ++i)
      {
    	const double xi = i*hh;
	file << std::setw(10) << std::setprecision(8) << xi << "   "
	     << std::setw(16) << std::setprecision(12) << alpha(xi)  << "   "
	     << std::setw(16) << std::setprecision(12) << beta(xi)   << "   "
	     << std::setw(16) << std::setprecision(12) << bdry_transfer(xi)  << "\n";
      }
    file.close ();

    std::string pc_file_name = "pc.txt";
    std::string full_pc_file_name = dir+"/pc.txt";

    std::ofstream pc_file (full_pc_file_name);
    pc_file << "#    S_w         pc_matrix(S_w)        pc_fracture(S_w)\n";
    for (unsigned int i = 0; i <= n; ++i)
    {
       	const double xi = i*hh;
       	if(xi < 0.1) continue;
   	    pc_file << std::setw(10) << std::setprecision(8) << xi << "   "
   	     << std::setw(16) << std::setprecision(12) << pc_matrix(xi)  << "   "
   	     << std::setw(16) << std::setprecision(12) << pc_fracture(xi)  << "\n";
    }
    pc_file.close ();

    std::string bdry_file_name = "bdry.txt";
    std::string full_bdry_file_name = dir+"/bdry.txt";

    std::ofstream bdry_file (full_bdry_file_name);
    bdry_file << "#    t          P(g(t))         g(t) \n";
    double dt = tend/n;
       for (unsigned int i = 0; i <= n; ++i)
       {
          	const double ti = i*dt;
      	    bdry_file << std::setw(10) << std::setprecision(8) << ti << "   "
      	     << std::setw(16) << std::setprecision(12) << bdry(ti) << "   "
			 << std::setw(16) << std::setprecision(12) << ptfun[function_index](ti)  << "\n";
       }
       bdry_file.close ();

       // Write gnuplot command file
    std::string file_name_gnu = dir + "/functions.gnu";
    std::ofstream out(file_name_gnu);
    auto tmp2 = aux::min_max(full_pc_file_name, 1);
    auto tmp3 = aux::min_max(full_pc_file_name, 2);
    auto max = std::max(tmp2.second, tmp3.second);
    out << "#set terminal epslatex size 3.5,2.62 color colortext\n";
    out << "#set terminal epslatex color # linewidth 2\n";
    out << "#set output \"fun-1-5-2.tex\"\n";

    out << "#set border linewidth 1\n";

    out << "set multiplot layout 2,2\n";

    out << "set xrange [0:1]\n";
    out << "set xtics ('$0$' 0.0, '$.2$' 0.2,'$.4$' 0.4,'$.6$' 0.6,'$.8$' 0.8, '$1$' 1.0) \n";

    out << "set xlabel '$S_w$'\n";
    out << "set ylabel '$P_c$ [bar]'\n";
    out << "set yrange [0.0:" << max << "]\n";
    out << "plot \"pc.txt\" u 1:2 w l title '$P_{c,m}(S_w)$',\\\n";
    out << "     \"pc.txt\" u 1:3 w l title '$P_{c,m}(S_w)$'\n";

    out << "unset ylabel\n";
    out << "set key right center\n";
    out << "#set grid\n";
    out << "set yrange [0.0:1]\n";
    out << "plot	 \"functions.txt\" u 1:4 w l title '${\\cal P}(S_w)$'\n";

    tmp2 = aux::min_max(full_bdry_file_name, 0);
    max = tmp2.second;

    out << "set ylabel \"time [days]\"\n";
    out << "set yrange [0:1.1]\n";
    out << "set key right bottom\n";
    out << "unset xlabel\n";
    out << "unset xrange\n";
    out << "unset xtics\n";
    out << "set xlabel 'time [days]'\n";
    out << "set xrange [0:" << max << "]\n";
    out << "set xtics\n";
    out << "plot \"bdry.txt\" u 1:2 w l title '${\\cal P}(S_w^f(t))$',\\\n";
    out << "     \"bdry.txt\" u 1:3 w l title '$S_w^f(t)$'\n";


    tmp2 = aux::min_max(full_file_name, 1);
    tmp3 = aux::min_max(full_file_name, 2);
    max = std::max(tmp2.second, tmp3.second);
    out << "unset ylabel\n";
    out << "unset xlabel\n";
    out << "unset xrange\n";
    out << "unset yrange\n";
    out << "set key right center\n";
    out << "set xlabel '$S_w$'\n";
    out << "set yrange [0:" << max << "]\n";
    out << "set xrange [0:1]\n";
    out << "set xtics ('$0$' 0.0, '$.2$' 0.2,'$.4$' 0.4,'$.6$' 0.6,'$.8$' 0.8, '$1$' 1.0) \n";
    out << "plot \"functions.txt\" u 1:2 w l title '$\\alpha_m(S_w)$',\\\n";
    out << "	 \"functions.txt\" u 1:3 w l title '$\\beta_m(S_w)$'\n";

    out << "unset multiplot\n";
    out << "pause -1\n";

	out.close();
}


void Params::gnu_compare_c(std::set<int> const & show_sim,
		                   std::string const & add_to_name) const{
	// total number of simulations to show
	int total_cnt = 0;
	for (unsigned int i = 0; i < size; ++i) if(simulation[i] and show_sim.count(i)) total_cnt++;
    if(total_cnt == 0) return;

//	const std::string & dir = params.date_and_time;
	std::string file = "flux"+add_to_name+".gnu";
	std::ofstream out(date_and_time + "/"+ file);
	out << "#set terminal epslatex size 3.5,2.62 color colortext\n";
	out << "#set terminal epslatex color # linewidth 2\n";
	out << "#set output \"flux-" << add_to_name << ".tex\"\n";

	out << "set xlabel \"time [days]\"\n";
	out << "set ylabel \"Nonwetting source [m^3/day]\"\n";
	out << "set key left center\n";
	out << "#set grid\n";

	std::pair<double, double> tmp={1.0E100, -1.0E100}; // min, max
	double time_max = 1.0;
	for (unsigned int i = 0; i < size; ++i) {
		if (simulation[i]) {
	     	std::string name = str_sname  + simulation_names[i] + "-flux.txt";

	     	auto tmp0 = aux::min_max(name, 0);
	     	time_max = tmp0.second;

	     	if(i == analytic_const || i == analytic_var ||
	     	   i == analytic_new ||i == analytic_new1){
	     	    auto tmp1 = aux::min_max(name, 1);
	     	    if(tmp1.first < tmp.first) tmp.first = tmp1.first;
	     	    if(tmp1.second > tmp.second) tmp.second = tmp1.second;
	     	}
	     	else{
	     	    auto tmp1 = aux::min_max(name, 1);
	     	    if(tmp1.first < tmp.first) tmp.first = tmp1.first;
	     	    if(tmp1.second > tmp.second) tmp.second = tmp1.second;
	     	    auto tmp2 = aux::min_max(name, 2);
	     	    if(tmp2.first < tmp.first) tmp.first = tmp2.first;
	     	    if(tmp2.second > tmp.second) tmp.second = tmp2.second;

	     	}
		}
	}

	out << "set xrange [0:" << time_max<< "]\n";

	double dy = tmp.second - tmp.first;
	assert(dy >= 0.0);
	tmp.first -= 0.05*dy;
	tmp.second += 0.05*dy;
	out << "set yrange [" << tmp.first <<":" << tmp.second << "]\n";
	out << "#set terminal postscript eps color solid lw 3\n";
	out << "#set output \"clin_Q.eps\"\n";
	int cnt = 0;
	for (unsigned int i = 0; i < size; ++i) {
     	std::string name = simulation_names[i];
		if (i == 0)
			out << "   plot \\\n";
		if (simulation[i] and show_sim.count(i)) {
			cnt++;
			if(i == analytic_const || i == analytic_var
					||i == analytic_new||i == analytic_new1)
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
			<< date_and_time << "; gnuplot " << file << "\n";
}
#endif /* SRC_PARAMETERS_HH_ */
