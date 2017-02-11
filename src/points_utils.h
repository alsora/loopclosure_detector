#pragma once
#include "defs.h"


namespace pr {

  //! appends to dest a sequence of points lying on a segment
  //! @param dest: the array od 3d points that will be augmented
  //! @param start: the first point of the segment
  //! @param end:   the last point on the segment
  //! @param density: how many points per meter
  void putPointsOnSegment3D(Vector3fVector& dest, 
			    const Eigen::Vector3f& start,
			    const Eigen::Vector3f& end,
			    float density);


  //! creates a random world of 3d points lying on a set of segments
  //! @param world_points: the array od 3d points that will be modified
  //! @param lower_left_bottom: one extrema of the bounding 
  //!                           box containing the world 
  //! @param upper_right_top:   other extrema of bounding box
  //! @param num_segments:      how many segments
  //! @param density:           how many points per meter in the lines
  void makeWorld(Vector3fVector& world_points,
		 const Eigen::Vector3f& lower_left_bottom,
		 const Eigen::Vector3f& upper_right_top,
		 int num_segments,
		 float density);



  //! puts in an int matrix the indices of 2d points.
  //! all cells in the matrix whose coordinates do not appear
  //! in the points will contain the invalid index -1.
  //! if points[i]=(r,c), the img[r,c]=i;
  //! if a point falls outside the image, it is dropped
  //! @param img: the image where to write (must be preallocated)
  //! @param points: the vector of points to render
  void putPointIndices(IntImage& img, 
		       const Vector2fVector& points);

  //! draws a set of 2d points on an image
  //! @param img: the (preallocated) dest image
  //! @param points: the array of points
  //! @param color: the color of the points
  //! @param radius: the size of the point in the image (in pixels)
  void drawPoints(RGBImage& img, 
		  const Vector2fVector& points, 
		  const cv::Scalar& color,
		  int radius, int thickness);
    
    void drawPoints(RGBImage& img,
                    const Vector2fVector& points,
                    const cv::Scalar& color,
                    int radius);



  //! draws a set of correspondences between 2d points
  //! @param img: the dest image (preallocated)
  //! @param reference_image_points: the first vector of points
  //! @param current_image_points: the second vector of points
  //! @param correspondences: the array of correspondences.
  //!        if correspondences[k] = (i,j), a line between
  //!        reference_image_points[i] and current_image_points[j]
  //!        will be drawn
  //! @param color: the color
  void drawCorrespondences(RGBImage& img,
			   const Vector2fVector& reference_image_points,
			   const Vector2fVector& current_image_points,
			   const IntPairVector& correspondences,
			   cv::Scalar color=cv::Scalar(0,255,0));


}
