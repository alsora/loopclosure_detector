#include "defs.h"
#include "database.h"
#include "utils.h"
#include "points_utils.h"
#include "solver.h"
#include "distance_map_correspondence_finder.h"
#include <unistd.h>



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;
using namespace pr;


int main (int argc, char** argv) {
    
    
    ifstream inputFile("../loop-detector-2d.txt");
    string odomFile("../trajectoryLaser.g2o");
    
    
    
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
    
    vector<vector<int>> groupedScans;
    vector<int> island;
    map<int, Eigen::Vector3f>::const_iterator it;
    
    for (int iD : sequence_list){
        it = posesMap.find(iD);
        if (it!= posesMap.end()){
            if (iD!=sequence_list[0]){
                groupedScans.push_back(island);  }
            island.clear();
            island.push_back(iD);
        }
        
        else {
            island.push_back(iD);
        }
    }
    
    
    
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
    float errorThreshold = 2.5;    //Maximum error between reference and target
    int correspondencesThreshold = 400;     //Minimum number of corresponding points between reference and target
    
    float islandThrehsold = 25;
    Solver solver;
    
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
    
    Eigen::Matrix3f transf;
    float chi;
    
    
    float currentIslandChi;
    Eigen::Matrix3f currentIslandTransf;
    int currentIslandFirstID;
    float bestChi;
    Eigen::Matrix3f bestTransf;
    int bestID;
    float bestIslandChi;
    Eigen::Matrix3f bestIslandTransf;
    int bestIslandFirstID;
    
    ofstream output_file;
    output_file.open("candidate_closures.txt");
    
    
    for (int referenceID : sequence_list){
    bestChi = errorThreshold + 1;
    float actualX = (posesMap.find(referenceID)->second)(0);
    cout<<referenceID<<endl;
        if (!searchClosure){    //I start looking for closures when the robot returns on its path
            
            searchClosure = (lastX > actualX);
            lastX = actualX;
            inversionScan = referenceID;
        }

           if (searchClosure){
    reference_points = cloudsMap.find(referenceID)->second;
    display_reference = database.pointsToDisplay(reference_points);
    ref_pose = posesMap.find(referenceID)->second;//Computation of current initial guess
    

    correspondence_finder.init(display_reference, //Init is necessary 1 time for each reference
                               rows,
                               cols,
                               max_distance);
    
    for (vector<int> island : groupedScans){
        bool first = true;
        for (int targetID : island){
            //int targetID = sequence_list[i];
            target_pose = posesMap.find(targetID)->second;
            float dista =abs(actualX - target_pose(0));
            
            if ((dista<threshold)&&(targetID < inversionScan)){
                
                iterate = true;
                
                
                
                
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
                
                new_points.clear();
                state = solver.getX();
                for (Eigen::Vector2f point : target_points){
                    
                    new_point = state.block<2,2>(0,0)*point + state.block<2,1>(0,2);
                    new_points.push_back(new_point);
                    
                }
                corrected_points = new_points; //Adjusted points according to current guess
                display_corrected = database.pointsToDisplay(corrected_points);
                correspondence_finder.compute(display_corrected);
                
                
                chi = solver.computeError(correspondence_finder.correspondences(),target_points,reference_points );
                
                
                if (first){
                    
                    first = false;
                    if ((chi >= islandThrehsold)||(correspondence_finder.correspondences().size()<correspondencesThreshold)){
                        break;  }
                    
                    else{
                        currentIslandChi = chi;
                        currentIslandTransf = state;
                        currentIslandFirstID = targetID;
                    }
                }
                
                if ((chi < bestChi)&&(correspondence_finder.correspondences().size()>=correspondencesThreshold)){
                    bestChi = chi;
                    bestTransf = state;
                    bestID = targetID;
                    bestIslandChi = currentIslandChi;
                    bestIslandTransf = currentIslandTransf;
                    bestIslandFirstID = currentIslandFirstID;
                }
                
                
                
                
            }
        }
    }
    
    
    if ((bestChi <= errorThreshold)){
        output_file<<"SCAN: "<<referenceID<<endl;
        output_file<<"BestIsland: "<< bestIslandFirstID << " with error " << bestIslandChi << " and transformation "<< mat2vec(bestIslandTransf).transpose() <<endl;
        output_file<<"BestCandidate: "<< bestID << " with error " << bestChi << " and transformation "<< mat2vec(bestTransf).transpose() <<endl;
    }
    }
    }
    
}
