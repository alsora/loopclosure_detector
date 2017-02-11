#pragma once
#include "defs.h"



namespace pr {
  class ScanDatabase{
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    ScanDatabase (int threshold);
          ScanDatabase ();
    
      int extractData (const std::string line);
      
      bool pointsFromScan (int sequence_number);
      
      Vector2fVector pointsToDisplay(Vector2fVector input_points);
      
      int extractOdometry(std::string file_name);
      
      std::map<int, std::vector<float>> getScanList ();
      
      std::map<int,Vector2fVector> getScanDrawings();
      
      std::map<int, Eigen::Vector3f> getOdometryList();

      
    
  protected:
      int _threshold = 30;
      bool first = true;
      std::string tag, topic, frame_id;
      float min_range, max_range, min_angle, max_angle, angle_increment;
      size_t scan_size;
      int n_rows = 500;
      int n_cols = 500;
      

      std::map<int, std::vector<float>> scan_list;
      std::map<int, Vector2fVector> scan_drawings;
      std::map<int, Eigen::Vector3f> odometry_list;
  };
}
