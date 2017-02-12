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
    
    
    //Conditions for being considered a candidate closure
    float errorThreshold = 2.5;    //Maximum error between reference and target
    int correspondencesThreshold = 400;     //Minimum number of corresponding points between reference and target
    float islandThrehsold = 25;     //Discard a whole island if its first element has an error bigger than this threshold
    int distanceThreshold = 5;      // Don't compare scans which are too far
    
    

    ScanDatabase database(distanceThreshold);
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
    
    //Get odometry from g2o file
    database.extractOdometry(odomFile);
    
    
    map<int,Vector2fVector> cloudsMap = database.getScanDrawings();
    map<int, Eigen::Vector3f> posesMap = database.getOdometryList();
    
    //Group scans into islands, according to their pose
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
    
    
    
    

    Solver solver;
    
    DistanceMapCorrespondenceFinder correspondence_finder;
    
    int rows = 400;     //Dimensions for the distance map
    int cols = 400;
    
    float max_distance=40;  //Threshold for nearest neighbors
    
    
    
    
    //Creation of all the variables needed during the execution of the program
    int inversionScanID;  //Stores the number of the scan after which we start looking for closures
    float lastX = 0;    //Stores the X coordinate of the previous pose
    bool searchClosure = false; //Until it's false, the system does not look for closures
    bool firstOfIsland; //Denotes that I'm comparing with the first scan of an island group
    bool iterate;               //Variable used inside the ICP iteration, to stop at convergence
    Eigen::Vector3f referencePose;   //Pose of the reference scan (the one for which I currently search closures)
    Eigen::Vector3f targetPose;   //Pose of the target scan (the current candidate for closure)
    Eigen::Matrix3f guess;      //initial guess stored as a matrix
    Eigen::Matrix3f state;        //Actual state of LS algorithm
    Vector2fVector referencePoints; //Vector containing the actual reference points
    Vector2fVector displayableReference; //referencePoints transformated to fit them into a matrix image
    Vector2fVector targetPoints;   //Vector containing the actual original target points
    Vector2fVector adjustedPoints;    //Vector containing the target points modified according actual state
    Vector2fVector displayableAdjusted;   //adjustedPoints transformated to fit them into a matrix image
    Vector2fVector newPoints;          //Vector used as storage for transforming points
    Eigen::Vector2f newPoint;          //Vector used as storage for coordinates
    
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
    output_file.open("../candidate_closures.txt");
    
    
    for (int referenceID : sequence_list){
        bestChi = errorThreshold + 1;
        float actualX = (posesMap.find(referenceID)->second)(0);
        cout<<referenceID<<endl;
        if (!searchClosure){    //I start looking for closures when the robot returns on its path
            
            searchClosure = (lastX > actualX);
            lastX = actualX;
            inversionScanID = referenceID;
        }

           if (searchClosure){
               referencePoints = cloudsMap.find(referenceID)->second;      //Get all information related to current reference scan
               displayableReference = database.pointsToDisplay(referencePoints);
               referencePose = posesMap.find(referenceID)->second;
    

               correspondence_finder.init(displayableReference, //Init is necessary 1 time for each reference scan
                               rows,
                               cols,
                               max_distance);
    
               for (vector<int> island : groupedScans){
                   firstOfIsland = true;
                   for (int targetID : island){
            
                       targetPose = posesMap.find(targetID)->second;
                       float distanceFromReference =abs(actualX - targetPose(0));
            
                       if ((distanceFromReference<distanceThreshold*1.5)&&(targetID < inversionScanID)){ // I run ICP only if the scan preceeds the reference one and it's not too far from it.
                
                           iterate = true;
                
                           guess = vec2mat(referencePose).inverse() * vec2mat(targetPose);
                
                
                           targetPoints = cloudsMap.find(targetID)->second;
                
                           solver.init(referencePoints, targetPoints,guess); //Initialization of the solver, 1 time for each target
                
                           while (iterate){
                    
                               newPoints.clear();
                               state = solver.getX();
                               for (Eigen::Vector2f point : targetPoints){
                                   newPoint = state.block<2,2>(0,0)*point + state.block<2,1>(0,2);
                                   newPoints.push_back(newPoint);
                               }
                               
                               adjustedPoints = newPoints; //Adjusted points according to current guess
                               displayableAdjusted = database.pointsToDisplay(adjustedPoints);
                               correspondence_finder.compute(displayableAdjusted);
                    
                               iterate = solver.oneRound(correspondence_finder.correspondences(),false);
                    
                           }
                
                           newPoints.clear();
                           state = solver.getX();
                           for (Eigen::Vector2f point : targetPoints){
                               newPoint = state.block<2,2>(0,0)*point + state.block<2,1>(0,2);
                               newPoints.push_back(newPoint);
                           }
                           
                           adjustedPoints = newPoints; //Adjusted points according to current guess
                           displayableAdjusted = database.pointsToDisplay(adjustedPoints);
                           correspondence_finder.compute(displayableAdjusted);
                
                
                           chi = solver.computeError(correspondence_finder.correspondences(),targetPoints,referencePoints );
                
                
                           if (firstOfIsland){//If the error of the first scan is too big, break goes directly to the next island. Otherwise I store the node information.
                    
                               firstOfIsland = false;
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
