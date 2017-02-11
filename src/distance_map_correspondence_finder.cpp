#include "distance_map_correspondence_finder.h"
#include "points_utils.h"

namespace pr {

  DistanceMapCorrespondenceFinder::DistanceMapCorrespondenceFinder(){
    _reference_points=0;
    _rows=0;
    _cols=0;
    _max_distance=0;
  }

  void DistanceMapCorrespondenceFinder::init(const Vector2fVector& reference_points, int rows_, int cols_, float max_distance){
    _max_distance=max_distance;
    _reference_points=&reference_points;
    _rows=rows_;
    _cols=cols_;
    
    _input_indices_image.create(_rows,_cols);
    _input_indices_image=-1;
    putPointIndices(_input_indices_image, *_reference_points);
    _distance_map.compute(_indices_image, 
			  _distance_image, 
			  _input_indices_image, 
			  _max_distance);
  }

  void DistanceMapCorrespondenceFinder::compute(const Vector2fVector& current_points) {
    _correspondences.resize(current_points.size());
    int num_correspondences=0;
    for (size_t current_idx=0; current_idx<current_points.size(); current_idx++){
      const Eigen::Vector2f current_point=current_points[current_idx];
      int r=current_point.y();
      int c=current_point.x();
      if (r<0||r>_rows)
	continue;
      if (c<0||c>_cols)
	continue;
      int reference_idx=_indices_image.at<int>(r,c);
      if (reference_idx<0)
	continue;
      IntPair& correspondence=_correspondences[num_correspondences];
      correspondence.first=reference_idx;
      correspondence.second=current_idx;
      num_correspondences++;
    }
    _correspondences.resize(num_correspondences);
  }

} // end namespace
