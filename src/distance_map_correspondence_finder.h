#pragma once
#include "distance_map.h"

namespace pr {

  /**
     Simple correspondence finder class that uses a distance map in 2D.
     To use it:
     - create an instance of the object
     - whenever the reference changes, call init. This will trigger a calculation
       of a new distance map
     - To check for associations, use the compute(...) function, passing the set of points
       for which you wish to find the corresponding reference
   */
  class DistanceMapCorrespondenceFinder{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    //! ctor
    DistanceMapCorrespondenceFinder();
    
    //! init method, constructs the distance map of size rows by cols, expanded at max_distance
    //! @param reference_points: array containing all 2d points of the model. Points outside the
    //! distance map are not rendered
    //! @param rows_: the how tall is the distance map
    //! @param cols_: how wide is the distance map
    //! @param max_distance: max expansion distance
    void init(const Vector2fVector& reference_points, int rows_, int cols_, float max_distance);
    
    //! computes the associations for the points passed as input. The result is internally stored
    //! @param current_points: the points for which you want to compute the associations
    void compute(const Vector2fVector& current_points);

    //! accessor to the last vector of correspondences computed by the finder
    inline const IntPairVector& correspondences() const {return _correspondences;}
    
    inline const FloatImage& distanceImage() const {return _distance_image;}
    inline float maxDistance() const {return _max_distance;}
  protected:
    float _max_distance;
    int _rows;
    int _cols;
    DistanceMap _distance_map;
    FloatImage _distance_image;
    IntImage _input_indices_image;
    IntImage _indices_image;
    const Vector2fVector* _reference_points;
    IntPairVector _correspondences;
  };

}//end namespace
