#include <Eigen/Cholesky>
#include <iostream>
#include "solver.h"
namespace pr {
  
  Solver::Solver(){
      
    _scan_points=0;
    _ref_points=0;
    _damping=1;
    _min_num_inliers=0;
    _num_inliers=0;
    _kernel_thereshold=2;       }
    

  void Solver::init(    const Vector2fVector& ref_points,
                        const Vector2fVector& scan_points,
                        const Eigen::Matrix3f initial_guess){
      
      _ref_points=&ref_points;
      _scan_points=&scan_points;
      X = initial_guess;
      oldX_1.setIdentity();
      oldX_2.setIdentity();
      _last_chi1 = 0;
      _last_chi2 = 0;
      _damping=1;
      _count = 0;       }
    
  
    
    float Solver::computeError (const IntPairVector& correspondences){
        
        float chi = 0;
        float chiTot = 0;
        float threshold = _kernel_thereshold;
        for (const IntPair& correspondence: correspondences){
            
            int ref_idx = correspondence.first;
            int obs_idx = correspondence.second;
        
            
            Eigen::Vector2f predicted_point = X.block<2,2>(0,0)*(*_scan_points)[obs_idx] + X.block<2,1>(0,2);
            
            Eigen::Vector2f error = predicted_point-(*_ref_points)[ref_idx];
            
            chi = error.dot(error);
            
            if (chi >threshold) {
                chiTot+=chi;    }
            
            else {
                chiTot+= chi;   }
                                        }
        
        return chiTot;
        
        
    }
    
    
    Eigen::Matrix3f Solver::computeH (const IntPairVector& correspondences){
        
        _H.setZero();
        
        for (const IntPair& correspondence: correspondences){
            int obs_idx = correspondence.second;
            
            Eigen::Vector2f predicted_point = X.block<2,2>(0,0)*(*_scan_points)[obs_idx] + X.block<2,1>(0,2);
            
            Matrix2_3f jacobian;
            jacobian.block<2,2>(0,0).setIdentity();
            
            Eigen::Vector2f vec(-predicted_point.y(),predicted_point.x());
            jacobian.block<2,1>(0,2)= vec;
            
            _H+=jacobian.transpose()*jacobian;
        }
        
        return _H;
        
    }
    
    
    

    
    void Solver::errorAndJacobian(Eigen::Vector2f& error,
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

    
    
  void Solver::linearize(const IntPairVector& correspondences, bool keep_outliers){
      _H.setZero();
      _b.setZero();
      _num_inliers=0;
      _chi_inliers=0;
      _chi_outliers=0;
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


          float lambda=1;
          bool is_inlier=true;
          if (chi>_kernel_thereshold){
              lambda=sqrt(_kernel_thereshold/chi);
              is_inlier=false;
              _chi_outliers+=chi;   }
          else {
              _chi_inliers+=chi;
              _num_inliers++;       }
      
          if (is_inlier || keep_outliers){
              _H+=J.transpose()*J*lambda;
              _b+=J.transpose()*e*lambda;   }
        }

  }

    
    
  bool Solver::oneRound(const IntPairVector& correspondences, bool keep_outliers){
      using namespace std;
      
      _count+=1;
      linearize(correspondences, keep_outliers);
      _H+=Eigen::Matrix3f::Identity()*_damping;
      
      if(_num_inliers<_min_num_inliers) {
          cerr << "too few inliers, skipping" << endl;
          return false; }
      
      //compute a solution
      Eigen::Vector3f dx = _H.ldlt().solve(-_b);
      X = v2t(dx)*X;
      
      //Compute error with the new solution
      float newChi = Solver::computeError (correspondences);

      if  (newChi < _chi_inliers){
          _damping *= 0.5;  }
      
      else (newChi >= _chi_inliers ){
          _damping *= 4;
          X = oldX_1;   }
      
      float convergenceThreshold = 0.05;
      
      if (_count > 10){
          convergenceThreshold = 0.2;
          }
      
      if ((abs(_last_chi1 - newChi)<convergenceThreshold)&&(abs(_last_chi2 - newChi)<convergenceThreshold)&&(abs(_last_chi1 - _last_chi2)<convergenceThreshold)){
          return false; } //Error has converged
      
      else {
          _last_chi2 = _last_chi1;
          _last_chi1 = newChi;
          oldX_2 = oldX_1;
          oldX_1 = X;
          return true;  }//Error has not converged, another round is necessary
  
    }
    

    
    Eigen::Matrix3f Solver::getX(){
        return X;   }
    
    
    Eigen::Matrix3f Solver::getH() {
        return _H;  }

}
