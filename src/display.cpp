#include "display.h"
#include <fstream>
namespace pr {
    ScanDatabase::ScanDatabase(int threshold){
        _threshold = threshold;
    }

    ScanDatabase::ScanDatabase(){
    }
    
    int ScanDatabase::extractData(const std::string line){
    std::string tag_n, topic_n, frame_id_n;
      int sequence_number_n;
      double timing_count_n;
      int odom1_n, odom2_n, odom3_n,odom4_n,odom5_n,odom6_n,odom7_n, odom8_n;
      float min_range_n, max_range_n, min_angle_n, max_angle_n, angle_increment_n, time_increment_n, scan_time_n;
      int scan_size_n;
      
    std::istringstream iss(line);
      
      
      iss >> tag_n >> topic_n >> frame_id_n>> sequence_number_n >> timing_count_n >>
      odom1_n >> odom2_n >> odom3_n >> odom4_n >> odom5_n >> odom6_n >> odom7_n >> odom8_n >> min_range_n >> max_range_n>> min_angle_n >> max_angle_n >> angle_increment_n >> time_increment_n >> scan_time_n >> scan_size_n;
      
      if (first){
          tag = tag_n;
          topic = topic_n;
          frame_id = frame_id_n;
          min_range = min_range_n;
          max_range = max_range_n;
          min_angle = min_angle_n;
          max_angle = max_angle_n;
          angle_increment = angle_increment_n;
          scan_size = scan_size_n;
          
          first = false;
          
      }
      
    std::vector<float> ranges(scan_size);
      
      
      for (int i=0; i<scan_size; i++){
          iss >> ranges[i];
      }
      
        scan_list.insert(std::pair<int, std::vector<float>>(sequence_number_n, ranges));

        return sequence_number_n;
  }

    
    
  
    bool ScanDatabase::pointsFromScan (int sequence_number){
        typedef std::map<int, std::vector<float>>::const_iterator it_type;
        it_type it = scan_list.find(sequence_number);
        if (it == scan_list.end()){
            std::cout<< sequence_number << " not found in scan_list"<< std::endl;
            return false;
        }
        else{
        int count = 0;
            Vector2fVector drawing;
        for (auto &val : it->second){
            if (val<_threshold){
            float angle = min_angle + count*angle_increment;
            float val1 =  val*cos(angle);
            float val2 = val*sin(angle);
            //std::cout <<val << " " << val1 << " " << val2 << std::endl;
            Eigen::Vector2f vec(val1,val2);
                drawing.push_back(vec);}
            count ++;
            
        }
            scan_drawings.insert(std::pair<int,Vector2fVector>(sequence_number, drawing));
            return true;
        }
        
    }
    
    

     
    Vector2fVector ScanDatabase:: pointsToDisplay(Vector2fVector input_points ){
        
        Vector2fVector corrected_points;
        
        for (Eigen::Vector2f point : input_points){
            
            float val1 =  n_cols/2 +(point(0))*max_range;
            float val2 =  n_rows/2 - point(1)*max_range ;
            
            Eigen::Vector2f vec(val1,val2);
            corrected_points.push_back(vec);
        }
        
        return corrected_points;
    }
    
    
    
    int ScanDatabase:: extractOdometry(std::string file_name){
        
        std::string line;
        std::string g2o_tag;
        Eigen::Vector3f pose;
        int correspondent_scan_tag = -1;
        float dummy;
        
        std::vector<float> ranges(scan_size);
        
        std::ifstream g2oFile(file_name);
        

        
        while (getline(g2oFile, line))
        {
            std::istringstream iss(line);
            iss >>g2o_tag;
            
            if (g2o_tag.compare("VERTEX_SE2")==0){
                iss>>correspondent_scan_tag >>pose(0) >> pose(1) >> pose(2);
            }
            else if (g2o_tag.compare("ROBOTLASER1")==0){
                iss >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy; // Neglect first 8 elements
                for (int i=0; i<scan_size; i++){
                    iss >> ranges[i];
                }
                typedef std::map<int,  std::vector<float>>::iterator it_type;
                for (it_type it = scan_list.begin(); it != scan_list.end(); ++it ){
                    if (it->second == ranges){
                        correspondent_scan_tag = it->first;
                        odometry_list.insert(std::pair<int,Eigen::Vector3f>(correspondent_scan_tag, pose));
                        
                        correspondent_scan_tag = -1;
                        pose.setZero();
                        
                        break; //Exit from the for loop
                    }
                }
            }
            else {
                std::cout<< "Invalid tag: "<<g2o_tag<<std::endl;
            }
            

        }

        
        return  odometry_list.size();
        
    }
    
    
    
    
    std::map<int,Vector2fVector> ScanDatabase::getScanDrawings(){
        return scan_drawings;
    }
    
    
    
    std::map<int, std::vector<float>> ScanDatabase::getScanList(){
        return scan_list;
    
    }
    
    std::map<int, Eigen::Vector3f> ScanDatabase::getOdometryList(){
        return odometry_list;
    }


}
