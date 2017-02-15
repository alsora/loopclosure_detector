#pragma once
#include <queue>
#include <boost/circular_buffer.hpp>
#include "defs.h"


namespace pr {

class Solver{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    //! ctor
    Solver();

    //! init method, call it at the beginning
    void init(
              const Vector2fVector& ref_points,
              const Vector2fVector& scan_points,
              const Eigen::Matrix3f initial_guess);
  
    inline float kernelThreshold() const {return _kernel_thereshold;}

    inline void setKernelThreshold(float kernel_threshold) 
    {_kernel_thereshold=kernel_threshold;}

    //! chi square of the "good" points
    const float chiInliers() const {return _chi_inliers;}
    
    //! chi square of the "bad" points
    const float chiOutliers() const {return _chi_outliers;}
    
    //! number of inliers (an inlier is a point whose error is below kernel threshold)
    const int numInliers() const {return _num_inliers;}
    
    //! performs one iteration of optimization
    //! param keep_outliers: if true, the outliers are considered in the optimization 
    //! (but cut by the kernel)
    bool oneRound(const IntPairVector& correspondences, bool keep_outliers);

    void predict(Eigen::Vector2f predicted_point, const Eigen::Vector2f correspondent_point);
      
    float computeError (const IntPairVector& correspondences);
    
    Eigen::Matrix3f computeH(const IntPairVector& correspondences);
    
    Eigen::Matrix3f getX();
    
    Eigen::Matrix3f getH();

    
  protected:
      



    void errorAndJacobian(Eigen::Vector2f& error,
                          Matrix2_3f& jacobian,
                          const Eigen::Vector2f& observed_point,
                          const Eigen::Vector2f& reference_point);

    void linearize(const IntPairVector& correspondences, bool keep_outliers);

    Eigen::Matrix3f X;              //state of the solver (SE(2))
    float _kernel_thereshold;        //< threshold for the kernel
    float _damping;                  //< damping, to slow the solution
    int _min_num_inliers;            //< if less inliers than this value, the solver stops
    const Vector2fVector* _ref_points;  //Storage for the reference points
    const Vector2fVector* _scan_points; //Storage for the observed points
    Eigen::Matrix3f _H;             //Storage for the H matrix
    Eigen::Vector3f _b;             //Storage for the b vector
    float _chi_inliers;
    float _chi_outliers;
    Eigen::Matrix3f oldX_1;         //Storage for the state at iteration i-1 (used for backup and for checking convergence)
    Eigen::Matrix3f oldX_2;         //Storage for the state at iteration i-2 (used for checking convergence)
    float _last_chi1;
    float _last_chi2;
    int _num_inliers;
    int _count;                     //Number of iterations already performed for the current clouds of points
      
  };

}
