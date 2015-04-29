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

/** MGF = Mesh generating function. Generate 1D mesh on interval [0,L]
 * adapted to a baoundary layer (BL) on left or both sides of [0,L].
 * Lecture Notes in Mathematics 1985
 * Torsten Linß
 *
 *  Layer-Adapted Meshes *  for Reaction-Convection-Diffusion Problems, 
 *  Springer
 *   /articles/Knjige/NumericalPDE
 *  [Torsten_Lin]_Layer-Adapted_Meshes_for_Reaction-C(BookZZ.org).pdf
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
      : L_(0.0), eps_(0.0), q_(0.0), sigma_(0.0), tau_(0.0),
        derivative_(0.0), iter_(0), N_(0)
      {
           L_     = params.L;
           eps_   = params.delta *params.delta;
           q_     = params.q;
           sigma_ = params.sigma;
           N_     = params.N;

           if(L_     <= 0.0)   throw std::runtime_error("L <= 0");
           if(eps_   <= 0.0)   throw std::runtime_error("eps <= 0");
           if(q_     <= 0.0)   throw std::runtime_error("q <= 0");
           if(sigma_ <= 0.0)   throw std::runtime_error("sigma <= 0");
           if(N_     <= 0  )   throw std::runtime_error("N <= 0");

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
          double yi =  t_to_x_general( i/static_cast<double>(N) );
          pts[i] = yi;
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
        // make N pair if it is already not.  
        if( (N/2)*2 != N) ++N;
        // Indices go as follows 
        // 0 1 2 ... M-1 M M+1 ...  N-1 N     where M = N/2. 
        int M = N/2;
        pts.resize(N+1);
        for(int i=0; i < M; ++i){
          double yi =  t_to_x_general( i/static_cast<double>(N) );
          double xi = yi; //(L_/2.0)*yi;
          pts[i] = xi;
          pts[N-i] = L_ - xi;
        }
        pts[M] = L_/2.0;
        return;
    }


  private:
    // Calculation of tau -- the end of BL part of the grid.
     void find_end_of_interval(){
        double tau_old = -100.0;
        double tau_next = 0.0;
        const double TOL = 1.0E-7;
        const int ITERMAX = 100;
        int iter = 0;
        while(std::abs(tau_old - tau_next) >= TOL and iter < ITERMAX)
        {
            tau_old = tau_next;
            tau_next = q_ -  (sigma_*eps_) *(0.5 - tau_old)/(L_/2 - t_to_x(tau_old));
            iter++;
//	    std::cout << tau_next << std::endl;
        }
        if(iter >= ITERMAX) 
          throw std::runtime_error("find_end_of_interval : maximim number of iterations excided");

        // save calculated values
        tau_  = tau_next;
        iter_ = iter;
        x_tau_ = t_to_x(tau_);
        derivative_ = der_t_to_x(tau_);
    }

     double x_to_t(double x) const { 
      return q_*(1-std::exp(-x/(sigma_*eps_)));
    }

     // transform parameter t \in [0,1] into grid point, but only in the
     // BL part of the grid
    double t_to_x(double t)  const { 
      if(t == 0.0) return 0.0;
      return -(sigma_*eps_)*std::log( (q_-t)/q_ );
    } 

    // Derivative of t_to_x(t)
    double der_t_to_x(double t)  const { 
      return (sigma_*eps_)*( 1.0/(q_-t) );
    } 

     // transform parameter t \in [0,1] into grid point, in both the
     // BL part of the grid, and the uniform part of the grid
    double t_to_x_general(double t) const {
      if(t <= tau_) return t_to_x(t);                          // curved part
      else          return x_tau_ + derivative_ * ( t - tau_); // uniform part
    }

    double L_;  // lenght of range interval
    double eps_;  
    double q_;   // q_ in (0,1)
    double sigma_; // grading of the mesh inside the layer
    double tau_;
    double x_tau_; 
    double derivative_;
    double iter_; // no of iterations to find tau
    int    N_;  //
};


#endif
