#include <iostream>
#include "points_utils.h"

namespace pr {
  using namespace std;


  void drawPoints(RGBImage& img, 
		  const Vector2fVector& points, 
		  const cv::Scalar& color,
		  int radius, int thickness){
    int rows=img.rows;
    int cols=img.cols;
    for (const Eigen::Vector2f point: points){

      int r=point.y();
      int c=point.x();
        
      if(r<0||r>rows)
	continue;
      if(c<0||c>cols)
	continue;
      cv::circle(img, cv::Point(c,r), radius, color,thickness);
    }
  } 

    void drawPoints(RGBImage& img,
                    const Vector2fVector& points,
                    const cv::Scalar& color,
                    int radius){
        int rows=img.rows;
        int cols=img.cols;
        for (const Eigen::Vector2f point: points){
            
            int r=point.y();
            int c=point.x();
            
            if(r<0||r>rows)
                continue;
            if(c<0||c>cols)
                continue;
            cv::circle(img, cv::Point(c,r), radius, color);
        }
    }
    
    
  void putPointIndices(IntImage& img, 
			const Vector2fVector& points){
      
    int rows=img.rows;
    int cols=img.cols;
    img=-1;
    for (size_t i=0; i<points.size(); i++){
      const Eigen::Vector2f& point=points[i];
      int r=point.y();
      int c=point.x();
      if(r<0||r>=rows)
	continue;
      if(c<0||c>=cols)
	continue;
      img.at<int>(r,c)=i;
    }
  }


  void drawCorrespondences(RGBImage& img,
			   const Vector2fVector& reference_image_points,
			   const Vector2fVector& current_image_points,
			   const IntPairVector& correspondences,
			   cv::Scalar color){
    for (const IntPair& correspondence: correspondences) {
      int ref_idx=correspondence.first;
      int curr_idx=correspondence.second;
      const Eigen::Vector2f& reference_point=reference_image_points[ref_idx];
      const Eigen::Vector2f& current_point=current_image_points[curr_idx];
      cv::line(img, 
	       cv::Point(reference_point.x(), reference_point.y()),
	       cv::Point(current_point.x(), current_point.y()),
	       color);
    }
  }


  void putPointsOnSegment3D(Vector3fVector& dest, 
		  const Eigen::Vector3f& start,
		  const Eigen::Vector3f& end,
		  float density){
    Eigen::Vector3f segment=end-start;
    size_t num_points=segment.norm()*density;
    if (num_points<1)
      num_points=1;
    segment *= (1.f/num_points);
    int k=dest.size();
    dest.resize(k+num_points);
    for (size_t i=0; i<num_points; i++, k++){
      dest[k]=start+segment*(float)i;
    }
  }

  void makeWorld(Vector3fVector& world_points,
		 const Eigen::Vector3f& lower_left_bottom,
		 const Eigen::Vector3f& upper_right_top,
		 int num_segments,
		 float density){
    for (int i=0; i<num_segments; i++){
      Eigen::Vector3f ranges=upper_right_top-lower_left_bottom;
      Eigen::Vector3f p0=
	lower_left_bottom+
	Eigen::Vector3f(ranges.x()*drand48(),
			ranges.y()*drand48(),
			ranges.z()*drand48());

      Eigen::Vector3f p1=
	lower_left_bottom+
	Eigen::Vector3f(ranges.x()*drand48(),
			ranges.y()*drand48(),
			ranges.z()*drand48());

      putPointsOnSegment3D(world_points, p0, p1, density);
    }
  }


}
