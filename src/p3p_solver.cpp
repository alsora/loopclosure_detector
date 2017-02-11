#include <Eigen/Cholesky>
#include <iostream>

#include "p3p_solver.h"
namespace pr {
  
  P3PSolver::P3PSolver(){
      
    _scan_points=0;
    _ref_points=0;
    _damping=1;
    _min_num_inliers=0;
    _num_inliers=0;
    _kernel_thereshold=10;
  }

  void P3PSolver::init(
	    const Vector2fVector& ref_points,
	    const Vector2fVector& scan_points,
                       const Eigen::Matrix3f initial_guess){
    _ref_points=&ref_points;
      _scan_points=&scan_points;
      X = initial_guess;
      oldX_1.setIdentity();
      oldX_2.setIdentity();
      _last_chi1 = 0;
      _last_chi2 = 0;
      
  }
  
    
    float P3PSolver::computeError (const IntPairVector& correspondences,
                                   const Vector2fVector observed_points,
                                   const Vector2fVector reference_points
                                   ){

        float chi = 0;
        int count = 0;
        float threshold = 10000;
        for (const IntPair& correspondence: correspondences){
            
            int ref_idx = correspondence.first;
            int obs_idx = correspondence.second;
        
            
             Eigen::Vector2f predicted_point = X.block<2,2>(0,0)*observed_points[obs_idx] + X.block<2,1>(0,2);
            
            Eigen::Vector2f error = predicted_point-reference_points[ref_idx];
            
            
            if (error.dot(error) >threshold) {
//std::cout<<"Outlier error: "<< error.dot(error) << " between point " <<(observed_points[obs_idx]).transpose()<<"("<< obs_idx << ") and " << (reference_points[ref_idx]).transpose()<<"("<<ref_idx <<")"<< std::endl;
                count++;

            }
            
            else {
                chi+= error.dot(error);
              //  std::cout<<error.dot(error)<<std::endl;
            }
     
            
        }
        //std::cout<<"Error: "<< chi<< " With "<< count<<" pairs with chi greater than "<< threshold << std::endl;
        return chi;
        
        
    }
    
    

    
    void P3PSolver::errorAndJacobian(Eigen::Vector2f& error,
                                     Matrix2_3f& jacobian,
                                     const Eigen::Vector2f& observed_point,
                                     const Eigen::Vector2f& reference_point){
        
        // compute the prediction
        Eigen::Vector2f predicted_point;
        predicted_point = X.block<2,2>(0,0)*observed_point + X.block<2,1>(0,2);
        // compute the error
        error=predicted_point-reference_point;
        // compute the jacobian of the transformation
        jacobian.block<2,2>(0,0).setIdentity();
        Eigen::Vector2f vec(-predicted_point.y(),predicted_point.x());
        jacobian.block<2,1>(0,2)= vec;
        
    }

  void P3PSolver::linearize(const IntPairVector& correspondences, bool keep_outliers){
    _H.setZero();
    _b.setZero();
    _num_inliers=0;
    _chi_inliers=0;
    _chi_outliers=0;
    //  std::cout<<"linearize"<<std::endl;
    for (const IntPair& correspondence: correspondences){
      Eigen::Vector2f e;
      Matrix2_3f J;
      int ref_idx=correspondence.first;
      int curr_idx=correspondence.second;
      errorAndJacobian(e,
				   J,
				   (*_scan_points)[curr_idx],
				   (*_ref_points)[ref_idx]);
        

      float chi=e.dot(e);
       // std::cout<<chi << " " ;

      float lambda=1;
      bool is_inlier=true;
      if (chi>_kernel_thereshold){
	lambda=sqrt(_kernel_thereshold/chi);
	is_inlier=false;
	_chi_outliers+=chi;
      } else {
	_chi_inliers+=chi;
	_num_inliers++;
      }
      
      if (is_inlier || keep_outliers){
	_H+=J.transpose()*J*lambda;
	_b+=J.transpose()*e*lambda;
      }
    }
              //std::cout<<std::endl;
  }

  bool P3PSolver::oneRound(const IntPairVector& correspondences, bool keep_outliers){
    using namespace std;
    linearize(correspondences, keep_outliers);
      _H+=Eigen::Matrix3f::Identity()*_damping;
    if(_num_inliers<_min_num_inliers) {
      cerr << "too few inliers, skipping" << endl;
      return false;
    }
    //compute a solution
      Eigen::Vector3f dx = _H.ldlt().solve(-_b);
      X = v2t(dx)*X;
      
      if ((_last_chi2!=0)&&(_last_chi2<_chi_inliers)){
          
          if (_last_chi1<_last_chi2){
       //  std::cout<<"STOP with "<< _last_chi1<<std::endl;
              X = oldX_1;
          }
          else {
           //   std::cout<<"STOP with "<< _last_chi2<<std::endl;
              X = oldX_2;
          }
          return false;
      }
      
      if ((std::abs(_last_chi1 - _chi_inliers)<0.15)&&(std::abs(_last_chi2 - _chi_inliers)<0.15)&&(std::abs(_last_chi1 - _last_chi2)<0.15)){
          
        //  std::cout<<"Convergence stop with "<< _chi_inliers<<std::endl;
          return false;
      }
      else {
          
          _last_chi2 = _last_chi1;
          _last_chi1 = _chi_inliers;
          oldX_2 = oldX_1;
          oldX_1 = X;
          return true;}
  }
    

    
    Eigen::Matrix3f P3PSolver::getX(){
        return X;
    }
}
