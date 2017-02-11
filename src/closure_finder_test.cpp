#include "defs.h"
#include "display.h"
#include "utils.h"
#include "points_utils.h"
#include "p3p_solver.h"
#include "distance_map_correspondence_finder.h"
#include <unistd.h>



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace pr;


int main (int argc, char** argv) {
    

    ifstream inputFile("../03-2DLoopDetector/loop-detector-2d.txt");
    string odomFile("../03-2DLoopDetector/trajectoryLaser.g2o");
    
    
   
    int threshold = 5;
    ScanDatabase database(threshold);
    vector<int> sequence_list;
    int scan_count = 0;
    
    string line;
    
    //Get image matrices from each scan line
    while (getline(inputFile, line))
    {
       int seq_number = database.extractData(line);
        sequence_list.push_back(seq_number);
        database.pointsFromScan(seq_number);
        scan_count++;
    }
    
    database.extractOdometry(odomFile);
    

   map<int,Vector2fVector> cloudsMap = database.getScanDrawings();
    map<int, Eigen::Vector3f> posesMap = database.getOdometryList();
    
    Eigen::Vector3f lastPose = posesMap.begin()->second; //initialize at first known pose
    
    for (int iD : sequence_list){
        if (posesMap.find(iD)!=posesMap.end()){ //element found
            lastPose = posesMap.find(iD)->second;
        }
        else {
            posesMap.insert(pair<int,Eigen::Vector3f>(iD, lastPose));
        }
    }



    
    //Conditions for being considered a candidate closure
    float errorThreshold = 0.015;    //Maximum error between reference and target
    int correspondencesThreshold = 400;     //Minimum number of corresponding points between reference and target
    
    
    P3PSolver solver;
    
    DistanceMapCorrespondenceFinder correspondence_finder;

    int rows = 400;     //Dimensions for the distance map
    int cols = 400;

    float max_distance=80;
    
    
    
    
    //Creation of all the variables needed during the execution of the program
    int inversionScan;  //Stores the number of the scan after which we start looking for closures
    float lastX = 0;    //Stores the X coordinate of the previous pose
    bool searchClosure = false; //Until it's false, the system does not look for closures
    bool iterate;               //Variable used inside the ICP iteration, to stop at convergence
    Eigen::Vector3f ref_pose;   //Pose of the reference scan (the one for which I currently search closures)
    Eigen::Vector3f target_pose;   //Pose of the target scan (the current candidate for closure)
    Eigen::Matrix3f guess;      //initial guess stored as a matrix
    Eigen::Matrix3f state;        //Actual state of LS algorithm
    Vector2fVector reference_points; //Vector containing the actual reference points
    Vector2fVector display_reference; //reference_points transformated to fit them into a matrix image
    Vector2fVector target_points;   //Vector containing the actual original target points
    Vector2fVector corrected_points;    //Vector containing the target points modified according actual state
    Vector2fVector display_corrected;   //corrected_points transformated to fit them into a matrix image
    Vector2fVector new_points;          //Vector used as storage for transforming points
    Eigen::Vector2f new_point;          //Vector used as storage for coordinates
    float best_chi;
    float best_target;
    Eigen::Vector3f best_transf;

    ofstream output_file;
    output_file.open("candidate_closures.txt");
    
    for (int j = 0; j<sequence_list.size();j++){ //Loop over all the scan (they will be considered as reference)
        best_chi = errorThreshold;
         int referenceID = sequence_list[j];
        float actualX = (posesMap.find(referenceID)->second)(0);
         cout<<referenceID<<endl;
        
        if (!searchClosure){    //I start looking for closures when the robot returns on its path
            
            searchClosure = (lastX > actualX);
            lastX = actualX;
            inversionScan = j;
        }
       reference_points = cloudsMap.find(referenceID)->second;
       display_reference = database.pointsToDisplay(reference_points);
        
        if (searchClosure){
            output_file<<"SCAN: "<< referenceID <<endl;
            correspondence_finder.init(display_reference, //Init is necessary 1 time for each reference
                                       rows,
                                       cols,
                                       max_distance);
            
            for (int i=0; i<inversionScan;i++){ //Loop over all the scans preceeding motion inversion
                int targetID = sequence_list[i];
                if (abs(actualX - (posesMap.find(targetID)->second)(0))<threshold){
                
            iterate = true;
                
                ref_pose = posesMap.find(referenceID)->second;//Computation of current initial guess
                target_pose = posesMap.find(targetID)->second;

                guess = vec2mat(ref_pose).inverse() * vec2mat(target_pose);
            
                               
                target_points = cloudsMap.find(targetID)->second;
                
                solver.init(reference_points, target_points,guess); //Initialization of the solver, 1 time for each target
                               
                while (iterate){
                    
                    new_points.clear();
                    state = solver.getX();
                    for (Eigen::Vector2f point : target_points){
                        
                        new_point = state.block<2,2>(0,0)*point + state.block<2,1>(0,2);
                        new_points.push_back(new_point);
                        
                    }
                    corrected_points = new_points; //Adjusted points according to current guess
                   display_corrected = database.pointsToDisplay(corrected_points);
                    correspondence_finder.compute(display_corrected);

                    iterate = solver.oneRound(correspondence_finder.correspondences(),false);
                
                }

                float chi = solver.computeError(correspondence_finder.correspondences(),target_points,reference_points );
                    
                  /*  if ((chi/correspondence_finder.correspondences().size() < best_chi)&&(correspondence_finder.correspondences().size()>=correspondencesThreshold)){
                        best_chi = chi/correspondence_finder.correspondences().size();
                        best_target = targetID;
                        best_transf = mat2vec(solver.getX());
                    }*/
                    
                if ((chi/correspondence_finder.correspondences().size()  <= errorThreshold)&&(correspondence_finder.correspondences().size()>=correspondencesThreshold)){
                    output_file<<"Candidate: "<< targetID << " with error " << chi/target_points.size() << " and transformation "<< mat2vec(solver.getX()).transpose() <<endl;

                }
                }
                
            }
            
         //   if (best_chi!=errorThreshold){
          //      output_file<<"SCAN "<< referenceID<< " BEST CANDIDATE: "<< best_target << " with mean error " << best_chi << " and transformation "<< best_transf.transpose() <<endl;
           // }
        }
        
        
        }

}
