#pragma once
#include <queue>
#include <boost/circular_buffer.hpp>
#include "defs.h"


namespace pr {

class P3PSolver{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    //! ctor
    P3PSolver();

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
    //! @param correspondences: the correspondences (first: measurement, second:model);
    //! param keep_outliers: if true, the outliers are considered in the optimization 
    //! (but cut by the kernel)
    bool oneRound(const IntPairVector& correspondences, bool keep_outliers);

         void predict(Eigen::Vector2f predicted_point, const Eigen::Vector2f correspondent_point);
      
      Eigen::Matrix3f getX();
      
      
      float computeError (const IntPairVector& correspondences,
                          const Vector2fVector reference_points,
                          const Vector2fVector observed_points);
      

      
  protected:
      



    void errorAndJacobian(Eigen::Vector2f& error,
			  Matrix2_3f& jacobian,
			  const Eigen::Vector2f& observed_point,
			  const Eigen::Vector2f& reference_point);

    void linearize(const IntPairVector& correspondences, bool keep_outliers);

      Eigen::Matrix3f X;
    float _kernel_thereshold;        //< threshold for the kernel
    float _damping;                  //< damping, to slow the solution
    int _min_num_inliers;            //< if less inliers than this value, the solver stops
    const Vector2fVector* _ref_points;
    const Vector2fVector* _scan_points;
      Eigen::Matrix3f _H;
      Eigen::Vector3f _b;
    float _chi_inliers;
    float _chi_outliers;
      Eigen::Matrix3f oldX_1;
      Eigen::Matrix3f oldX_2;
      float _last_chi1;
      float _last_chi2;
    int _num_inliers;
      
  };

}
