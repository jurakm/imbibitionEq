#ifndef _MGF_HH_INCLUDED__
#define _MGF_HH_INCLUDED__

// Mesh generation function
//
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

/** \brief Mesh generating function.
 *
 * Generate 1D mesh on interval [0,L]
 * adapted to a boundary layer (BL) on left side or on the both sides of [0,L].
 * Lecture Notes in Mathematics 1985
 * Torsten Linß
 *
 *  Layer-Adapted Meshes *  for Reaction-Convection-Diffusion Problems, 
 *  Springer
 *   /articles/Knjige/NumericalPDE
 *  [Torsten_Lin]_Layer-Adapted_Meshes_for_Reaction-C(BookZZ.org).pdf
 *
 *  TODO In the case of one sided boundary layer the grid is not smooth at x = L/2.
 * */
template <typename Params>
class MGF{
  public:
    /**
     * @param L    = length of the domain
     * @param eps  = small parameter in front of second derivative in PDE
     * @param q   = portion of the mesh used to resolve the layer  (0< q < 1)
     * @param sigma = grading of the mesh inside the layer (sigma > 0)
     * @param beta = slope of the boundary layer (>0, an estimate)
     **/
    MGF(Params const & params)
      : L_(0.0),  q_(0.0), C_(0.0), tau_(0.0), x_tau_(0.0), derivative_(0.0), iter_(0), N_(0)
      {
           L_     = params.L - params.delta;  // 
           double eps   = params.scaled_delta; // *params.delta;
           q_     = params.q;
           double sigma = params.sigma;
           N_     = params.N;
           C_     = eps * sigma;

           if(L_ <= 0.0)   throw std::runtime_error("L <= 0");
           if(C_ <= 0.0)   throw std::runtime_error("eps*sigma <= 0");
           if(q_ <= 0.0)   throw std::runtime_error("q <= 0");
           if(q_ >= 0.5)   throw std::runtime_error("q >= 0.5! Keep q in (0,0.5).");
           if(N_ <= 0  )   throw std::runtime_error("N <= 0");

           if(C_/q_ >= L_)
        	   std::cout << "MGF: warning. Boundary layer is not detected.\n";
           find_end_of_interval();
      }

    /** On [0,tau] generating function is curved. */   
    double tau() const  { return tau_; }
    /** Number of iterations to calculate tau.   */
    int    iter() const { return iter_;}

    /**  Cover interval [0,L] with a grid of N+1 points 
     *   assuming that the boundary layer is on the left side of the interval.
     *   @param N = number of points in the grid
     *   @param pts = vector of calculated points.
     */
    void one_side_interval(std::vector<double> & pts){
        int N = N_;
        if( N <= 0 )   throw std::string("N <= 0");
        pts.resize(N+1);
        for(int i=0; i <= N; ++i){
        	const double t_i = i/static_cast<double>(N); // t_i \in [0,1]
            const double y_i =  t_to_x_general( t_i );
            pts[i] = y_i;
        }
        return;
    }
  
    /**  Cover interval [0,L] with a grid of N+1 (or N+2 if N is impair) points 
     *   assuming that the boundary layer is on the both sides of the interval.
     *   @param N = number of points in the grid
     *   @param pts = vector of calculated points.
     */
    void double_side_interval(std::vector<double> & pts){
        int N = N_;
        if( N <= 0 )   throw std::string("N <= 0");
        // make N pair if it is not already.
        if( (N/2)*2 != N) ++N;
        // Indices go as follows 
        // 0 1 2 ... M-1 M M+1 ...  N-1 N     where M = N/2. 
        int M = N/2;
        pts.resize(N+1);
        for(int i=0; i < M; ++i){
          double t_i = i/static_cast<double>(N); // t_i \in [0,1/2)
          double x_i =  t_to_x_general( t_i  );
          pts[i] = x_i;
          pts[N-i] = L_ - x_i;
        }
        pts[M] = L_/2.0;
        return;
    }

  private:
    // Calculation of tau -- the end of BL part of the grid.
     void find_end_of_interval(){
    	if(C_/q_ >= L_){  // there is no boundary layer
    		tau_ = 0.0;
    		x_tau_ = 0.0;
    		derivative_ = L_;
    	}
        double tau_old = -100.0;
        double tau_next = 0.0;
        const double TOL = 1.0E-7;
        const int ITERMAX = 100;
        int iter = 0;
        while(std::abs(tau_old - tau_next) >= TOL and iter < ITERMAX)
        {
            tau_old = tau_next;
            tau_next = q_ -  C_ *(0.5 - tau_old)/(L_/2 - t_to_x(tau_old));
            iter++;
//	    std::cout << tau_next << std::endl;
        }
        if(iter >= ITERMAX) 
          throw std::runtime_error("find_end_of_interval : maximum number of iterations exceeded");

        // save calculated values
        tau_  = tau_next;
        iter_ = iter;
        x_tau_ = t_to_x(tau_);
        derivative_ = der_t_to_x(tau_);

//        const double residual = derivative_ - (L_/2 - x_tau_)/(0.5-tau_);
//        std::cout << "MGF::find_end_of_interval: residual = " << residual << std::endl;
    }

     double x_to_t(double x) const { 
      return q_*(1-std::exp(-x/(C_)));
    }

     // transform parameter t \in [0,1] into grid point, but only in the
     // BL part of the grid
    double t_to_x(double t)  const {
      assert( t < q_ );
      if(t == 0.0) return 0.0;
      return -C_*std::log( (q_-t)/q_ );
    } 

    // Derivative of t_to_x(t)
    double der_t_to_x(double t)  const { 
      return C_*( 1.0/(q_-t) );
    } 

     // transform parameter t \in [0,1] into grid point, in both the
     // BL part of the grid, and the uniform part of the grid
    double t_to_x_general(double t) const {
      if(t <= tau_) return t_to_x(t);                          // curved part
      else          return x_tau_ + derivative_ * ( t - tau_); // uniform part
    }

    double L_;  // length of range interval
//    double eps_;
    double q_;   // q_ in (0,1)
//    double sigma_; // grading of the mesh inside the layer
    double C_;
    double tau_;
    double x_tau_; 
    double derivative_;
    double iter_; // no of iterations to find tau
    int    N_;  //
};


#endif
